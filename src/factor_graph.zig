const std = @import("std");
const blasfeo = @import("blasfeo.zig");

pub fn removeConst(comptime T: type, ptr: *const T) *T {
    return @intToPtr(*T, @ptrToInt(ptr));
}

pub fn allocVecs(allocator: std.mem.Allocator, dim: c_int, dvecs: []*blasfeo.Dvec) !void {
    for (dvecs) |dvec| {
        dvec.mem = 0;
    }
    errdefer {
        for (dvecs) |dvec| {
            if (dvec.mem != 0) {
                allocator.destroy(dvec.mem);
            }
        }
    }
    const memsize = blasfeo.memsize_vec(dim);
    for (dvecs) |dvec| {
        var bytes = try allocator.alignedAlloc(u8, blasfeo.CACHE_LINE_SIZE, memsize);
        blasfeo.create_dvec(dim, dvec, bytes.ptr);
    }
}

pub fn allocMats(allocator: std.mem.Allocator, dim: c_int, dmats: []*blasfeo.Dmat) !void {
    for (dmats) |dmat| {
        dmat.mem = 0;
    }
    errdefer {
        for (dmats) |dmat| {
            if (dmat.mem != 0) {
                allocator.destroy(dmat.mem);
            }
        }
    }
    const memsize = blasfeo.memsize_mat(dim, dim);
    for (dmats) |dmat| {
        var bytes = try allocator.alignedAlloc(u8, blasfeo.CACHE_LINE_SIZE, memsize);
        blasfeo.create_dmat(dim, dim, dmat, bytes.ptr);
    }
}

var ARENA = blasfeo.Arena.init(std.heap.c_allocator);

pub fn infoToCov(info_vec: blasfeo.Vector, info_mat: blasfeo.Matrix, mean: blasfeo.Vector, cov: blasfeo.Matrix) !void {
    try ARENA.pushEnv();
    defer ARENA.popEnv();

    var chol = try ARENA.matrix(cov.rows, cov.rows);
    blasfeo.potrf_l(info_mat, chol);

    blasfeo.gese(0.0, cov);
    blasfeo.diare(1.0, cov);
    blasfeo.trsm_llnn(1.0, chol, cov, cov);
    blasfeo.trsm_lltn(1.0, chol, cov, cov);

    blasfeo.trsv_lnn(chol, info_vec, mean);
    blasfeo.trsv_ltn(chol, mean, mean);
}

pub fn covToInfo(mean: blasfeo.Vector, cov: blasfeo.Matrix, info_vec: blasfeo.Vector, info_mat: blasfeo.Matrix) !void {
    // actually this is the same function as infoToCov
    try infoToCov(mean, cov, info_vec, info_mat);
}

const FactorGraph = struct {
    const Self = @This();

    // zig fmt: off
    const Variable = struct {
        id: usize,
        dim : c_int,
        cov_det : f64 = 0,
        kl_divergence : f64 = 0,
        mean_dvec: blasfeo.Dvec = undefined,
        cov_dmat: blasfeo.Dmat = undefined,
        info_dvec: blasfeo.Dvec = undefined,
        info_dmat: blasfeo.Dmat = undefined,
        prior_info_dvec: blasfeo.Dvec = undefined,
        prior_info_dmat: blasfeo.Dmat = undefined,
        first_factor_slot : ?usize = null,
        last_factor_slot : ?usize = null,

        pub fn mean(self: *Variable) blasfeo.Vector {
            return blasfeo.Vector.attach(&self.mean_dvec);
        }
        pub fn cov(self: *Variable) blasfeo.Matrix {
            return blasfeo.Matrix.attach(&self.cov_dmat);
        }
        pub fn info_vec(self: *Variable) blasfeo.Vector {
            return blasfeo.Vector.attach(&self.info_dvec);
        }
        pub fn info_mat(self: *Variable) blasfeo.Matrix {
            return blasfeo.Matrix.attach(&self.info_dmat);
        }
        pub fn prior_info_vec(self: *Variable) blasfeo.Vector {
            return blasfeo.Vector.attach(&self.prior_info_dvec);
        }
        pub fn prior_info_mat(self: *Variable) blasfeo.Matrix {
            return blasfeo.Matrix.attach(&self.prior_info_dmat);
        }

        pub fn init(self: *Variable, prior_mean : []const f64, prior_cov : []const f64) !void {
            const dim = self.dim;
            std.debug.assert(prior_mean.len == dim);
            std.debug.assert(prior_cov.len == dim*dim);
            
            try ARENA.pushEnv();
            defer ARENA.popEnv();

            // get matrix and vector views into variable data
            var vmean = self.mean();
            var vcov = self.cov();
            var vinfo_vec = self.info_vec();
            var vinfo_mat = self.info_mat();
            var vprior_info_vec = self.prior_info_vec();
            var vprior_info_mat = self.prior_info_mat();

            // write mean and cov
            vcov.setFromSlice(prior_cov);
            vmean.setFromSlice(prior_mean);

            // set prior info vec/mat
            try covToInfo(vmean, vcov, vprior_info_vec, vprior_info_mat);

            // copy values into the current info and prior vec
            blasfeo.veccp(vprior_info_vec, vinfo_vec);
            blasfeo.gecp(vprior_info_mat, vinfo_mat);
        }
    };

    const Factor = struct {
        id: usize,
        dim: c_int,
        info_dvec : blasfeo.Dvec = undefined,
        info_dmat : blasfeo.Dmat = undefined,
        info_base_dvec : blasfeo.Dvec = undefined,
        info_base_dmat : blasfeo.Dmat = undefined,
        
        slot_begin : usize = 0,
        slot_end : usize = 0,

        pub fn info_mat(self: *Factor) blasfeo.Matrix {
            return blasfeo.Matrix.attach(&self.info_dmat);
        }

        pub fn info_vec(self: *Factor) blasfeo.Vector {
            return blasfeo.Vector.attach(&self.info_dvec);
        }
        
        pub fn info_base_mat(self: *Factor) blasfeo.Matrix {
            return blasfeo.Matrix.attach(&self.info_base_dmat);
        }

        pub fn info_base_vec(self: *Factor) blasfeo.Vector {
            return blasfeo.Vector.attach(&self.info_base_dvec);
        }

        pub fn finishInit(self: *Factor) void {
            blasfeo.gecp(self.info_base_mat(), self.info_mat());
            blasfeo.veccp(self.info_base_vec(), self.info_vec());
        }

        pub fn initWithSlice(self: *Factor, info_vec_data : []const f64, info_mat_data : []const f64) void {
            self.info_base_mat().setFromSlice(info_mat_data);
            self.info_base_vec().setFromSlice(info_vec_data);
            self.finishInit();
        }

        pub fn init(self: *Factor, info_vec_data: blasfeo.Vector, info_mat_data : blasfeo.Matrix) void {
            blasfeo.gecp(info_mat_data, self.info_base_mat());
            blasfeo.veccp(info_vec_data, self.info_base_vec());
            self.finishInit();
        }

        pub fn beginRelinearize(self: *Factor) void {
            blasfeo.gead(-1.0, self.info_base_mat(), self.info_mat());
            blasfeo.vecad(-1.0, self.info_base_vec(), self.info_vec());
        }

        pub fn endRelinearize(self: *Factor) void {
            blasfeo.gead(1.0, self.info_base_mat(), self.info_mat());
            blasfeo.vecad(1.0, self.info_base_vec(), self.info_vec());
        }

        pub fn relinearizeWithSlice(self: *Factor, info_vec_data : []const f64, info_mat_data : []const f64) void {
            self.beginRelinearize();
            self.info_base_mat().setFromSlice(info_mat_data);
            self.info_base_vec().setFromSlice(info_vec_data);
            self.endRelinearize();
        }

        pub fn relinearize(self: *Factor, info_vec_data : blasfeo.Vector, info_mat_data : blasfeo.Matrix) void {
            self.beginRelinearize();
            blasfeo.gecp(info_mat_data, self.info_base_mat());
            blasfeo.veccp(info_vec_data, self.info_base_vec());
            self.endRelinearize();
        }        
    };

    const FactorSlot = struct {
        const Child = @This();
        
        factor_id : usize,
        variable_id : usize,
        dim : c_int,
        dim_offset : c_int,
        next_same_variable_slot : ?usize = null,
        to_variable_dvec : blasfeo.Dvec = undefined,
        to_variable_dmat : blasfeo.Dmat = undefined,
        from_variable_dvec : blasfeo.Dvec = undefined,
        from_variable_dmat : blasfeo.Dmat = undefined,

        pub fn to_variable_vec(self: *Child) blasfeo.Vector {
            return blasfeo.Vector.attach(&self.to_variable_dvec);
        }

        pub fn to_variable_mat(self: *Child) blasfeo.Matrix {
            return blasfeo.Matrix.attach(&self.to_variable_dmat);
        }

        pub fn from_variable_vec(self: *Child) blasfeo.Vector {
            return blasfeo.Vector.attach(&self.from_variable_dvec);
        }

        pub fn from_variable_mat(self: *Child) blasfeo.Matrix {
            return blasfeo.Matrix.attach(&self.from_variable_dmat);
        }
    };

    const VariableFactorIterator = struct {
        const Parent = FactorGraph;
        const Child = @This();

        parent : *FactorGraph,
        variable: *Variable,
        current_slot_idx : ?usize = null,
        done: bool = false,

        pub fn next(self: *Child) ?*Factor {
            if (self.done) return null;

            // advance or initialize
            if (self.current_slot_idx) |slot_idx| {
                // advance
                self.current_slot_idx = self.parent.factor_slots.items[slot_idx].next_same_variable_slot;
            } else if (self.variable.first_factor_slot) |first_slot_idx| {
                // initialize
                self.current_slot_idx = first_slot_idx;
            }

            // return corresponding factor
            if (self.current_slot_idx) |slot_idx| {
                const factor_id = self.parent.factor_slots.items[slot_idx].factor_id;
                return &self.parent.factors.items[factor_id];
            } else {
                self.done = true;
                return null;
            }
        }
    };

    // zig fmt: on

    allocator: std.mem.Allocator,
    factors: std.ArrayList(Factor),
    factor_slots: std.ArrayList(FactorSlot),
    variables: std.ArrayList(Variable),

    pub fn init(allocator: std.mem.Allocator) Self {
        var self: Self = undefined;
        self.allocator = allocator;
        self.factors = std.ArrayList(Factor).init(allocator);
        self.factor_slots = std.ArrayList(FactorSlot).init(allocator);
        self.variables = std.ArrayList(Variable).init(allocator);
        return self;
    }

    pub fn deinit(self: *Self) void {
        self.factors.deinit();
        self.factor_slots.deinit();
        self.variables.deinit();
    }

    pub fn variableFactorsIterator(self: *Self, variable_id: usize) VariableFactorIterator {
        var it = VariableFactorIterator{ .variable = &self.variables.items[variable_id], .parent = self };
        return it;
    }

    pub fn addVariable(self: *Self, dim: c_int) !*Variable {
        var variable = Variable{ .id = self.variables.items.len, .dim = dim };
        var vecs: [3]*blasfeo.Dvec = .{ &variable.mean_dvec, &variable.info_dvec, &variable.prior_info_dvec };
        try allocVecs(self.allocator, dim, &vecs);
        errdefer {
            for (vecs) |vec| {
                self.allocator.destroy(vec.mem);
            }
        }
        var mats: [3]*blasfeo.Dmat = .{ &variable.cov_dmat, &variable.info_dmat, &variable.prior_info_dmat };
        try allocMats(self.allocator, dim, &mats);
        errdefer {
            for (mats) |mat| {
                self.allocator.destroy(mat.mem);
            }
        }

        try self.variables.append(variable);
        return &self.variables.items[self.variables.items.len - 1];
    }

    pub fn addFactor(self: *Self, variable_ids: []usize) !*Factor {
        // std.debug.print("Adding factor with ids {d}\n", .{variable_ids});
        // factor dim is the sum of the variable dims
        var dim: c_int = 0;
        for (variable_ids) |variable_id| {
            dim += self.variables.items[variable_id].dim;
        }
        var factor = Factor{ .id = self.factors.items.len, .dim = dim };

        var factor_vecs: [2]*blasfeo.Dvec = .{ &factor.info_dvec, &factor.info_base_dvec };
        try allocVecs(self.allocator, dim, &factor_vecs);
        errdefer {
            for (factor_vecs) |factor_vec| {
                self.allocator.destroy(factor_vec.mem);
            }
        }
        blasfeo.vecse(0.0, factor.info_vec());
        blasfeo.vecse(0.0, factor.info_base_vec());

        var factor_mats: [2]*blasfeo.Dmat = .{ &factor.info_dmat, &factor.info_base_dmat };
        try allocMats(self.allocator, dim, &factor_mats);
        errdefer {
            for (factor_mats) |factor_mat| {
                self.allocator.destroy(factor_mat.mem);
            }
        }
        blasfeo.gese(0.0, factor.info_mat());
        blasfeo.gese(0.0, factor.info_base_mat());

        const slot_begin = self.factor_slots.items.len;
        errdefer {
            while (self.factor_slots.items.len > slot_begin) {
                var slot = self.factor_slots.pop();
                self.allocator.destroy(slot.to_variable_dvec.mem);
                self.allocator.destroy(slot.from_variable_dvec.mem);
                self.allocator.destroy(slot.to_variable_dmat.mem);
                self.allocator.destroy(slot.from_variable_dmat.mem);
            }
        }

        var dim_offset: c_int = 0;
        for (variable_ids) |variable_id| {
            const variable = &self.variables.items[variable_id];
            var slot = FactorSlot{
                .factor_id = factor.id,
                .variable_id = variable_id,
                .dim = variable.dim,
                .dim_offset = dim_offset,
            };
            var vecs: [2]*blasfeo.Dvec = .{ &slot.to_variable_dvec, &slot.from_variable_dvec };
            try allocVecs(self.allocator, variable.dim, &vecs);
            errdefer {
                for (vecs) |vec| {
                    self.allocator.destroy(vec.mem);
                }
            }
            var mats: [2]*blasfeo.Dmat = .{ &slot.to_variable_dmat, &slot.from_variable_dmat };
            try allocMats(self.allocator, variable.dim, &mats);
            errdefer {
                for (mats) |mat| {
                    self.allocator.destroy(mat.mem);
                }
            }

            // send initial message to this factor / slot
            blasfeo.gecp(variable.info_mat(), factor.info_mat().block(dim_offset, dim_offset, variable.dim, variable.dim));
            blasfeo.veccp(variable.info_vec(), factor.info_vec().segment(dim_offset, variable.dim));
            blasfeo.gecp(variable.info_mat(), slot.from_variable_mat());
            blasfeo.veccp(variable.info_vec(), slot.from_variable_vec());

            try self.factor_slots.append(slot);
            dim_offset += variable.dim;
        }
        const slot_end = self.factor_slots.items.len;

        factor.slot_begin = slot_begin;
        factor.slot_end = slot_end;

        try self.factors.append(factor);

        // update internal linked lists inside the factor_slots
        // note:
        //   this is easier to do after we know no more errors can occur
        //   since unwinding is harder. do not try to merge this with the
        //   previous loop
        var slot_idx = slot_begin;
        while (slot_idx < slot_end) : (slot_idx += 1) {
            var slot = &self.factor_slots.items[slot_idx];
            var variable = &self.variables.items[slot.variable_id];
            if (variable.last_factor_slot) |prior_slot_idx| {
                // this variable is used in a prior slot
                var prior_slot = &self.factor_slots.items[prior_slot_idx];
                prior_slot.next_same_variable_slot = slot_idx;
                variable.last_factor_slot = slot_idx;
            } else {
                // this variable is not used in a prior slot
                variable.first_factor_slot = slot_idx;
                variable.last_factor_slot = slot_idx;
            }
        }

        return &self.factors.items[factor.id];
    }

    pub fn propagateVariable(self: *Self, variable_id: usize) !void {
        const variable = &self.variables.items[variable_id];
        if (variable.first_factor_slot == null) {
            // this variable is not involved with any factors
            return;
        }

        // std.debug.print("Propagating var {d}\n", .{variable.id});

        var next_slot: ?*FactorSlot = &self.factor_slots.items[variable.first_factor_slot.?];
        while (next_slot != null) {
            var slot = next_slot.?;
            if (slot.next_same_variable_slot) |next_slot_idx| {
                next_slot = &self.factor_slots.items[next_slot_idx];
            } else {
                next_slot = null;
            }

            try ARENA.pushEnv();
            defer ARENA.popEnv();

            var factor = &self.factors.items[slot.factor_id];
            // std.debug.print("Sending message to factor {d}\n", .{factor.id});

            // update the to-factor message from this variable
            // this is the variable info minus the info that it got from this slot
            blasfeo.gecp(variable.info_mat(), slot.from_variable_mat());
            blasfeo.gead(-1.0, slot.to_variable_mat(), slot.from_variable_mat());
            blasfeo.veccp(variable.info_vec(), slot.from_variable_vec());
            blasfeo.vecad(-1.0, slot.to_variable_vec(), slot.from_variable_vec());

            // std.debug.print("Info vec\n", .{});
            // blasfeo.print_vec(slot.from_variable_vec());
            // std.debug.print("Info mat\n", .{});
            // blasfeo.print_mat(slot.from_variable_mat());

            // incorporate the to-factor message into the factor info
            var info_block = factor.info_mat().block(slot.dim_offset, slot.dim_offset, variable.dim, variable.dim);
            var info_base_block = factor.info_base_mat().block(slot.dim_offset, slot.dim_offset, variable.dim, variable.dim);
            blasfeo.gecp(info_base_block, info_block);
            blasfeo.gead(1.0, slot.from_variable_mat(), info_block);

            var info_segment = factor.info_vec().segment(slot.dim_offset, variable.dim);
            var info_base_segment = factor.info_base_vec().segment(slot.dim_offset, variable.dim);
            blasfeo.veccp(info_base_segment, info_segment);
            blasfeo.vecad(1.0, slot.from_variable_vec(), info_segment);

            // std.debug.print("Updated Factor Info vec\n", .{});
            // blasfeo.print_vec(factor.info_vec());
            // std.debug.print("Updated Factor Info mat\n", .{});
            // blasfeo.print_mat(factor.info_mat());

            // marginalize out every other slot
            var work_mat = try ARENA.matrix(factor.dim, factor.dim);
            var work_mat_inv = try ARENA.matrix(factor.dim, factor.dim);
            var work_vec = try ARENA.vector(factor.dim);

            var other_slot_idx = factor.slot_begin;
            while (other_slot_idx < factor.slot_end) : (other_slot_idx += 1) {
                var other_slot = &self.factor_slots.items[other_slot_idx];
                var other_dim = other_slot.dim;

                blasfeo.gecp(factor.info_mat(), work_mat);
                blasfeo.veccp(factor.info_vec(), work_vec);

                // remove contribution from other slot
                // zig fmt: off
                blasfeo.gecp(
                    factor.info_base_mat().block(other_slot.dim_offset, other_slot.dim_offset, other_dim, other_dim),
                    work_mat.block(other_slot.dim_offset, other_slot.dim_offset, other_dim, other_dim)
                );
                blasfeo.veccp(
                    factor.info_base_vec().segment(other_slot.dim_offset, other_dim),
                    work_vec.segment(other_slot.dim_offset, other_dim)
                );

                try infoToCov(work_vec, work_mat, work_vec, work_mat_inv);
                
                // retreive the block
                var work_block = work_mat_inv.block(other_slot.dim_offset, other_slot.dim_offset, other_dim, other_dim);
                var work_segment = work_vec.segment(other_slot.dim_offset, other_dim);
                try covToInfo(work_segment, work_block, other_slot.to_variable_vec(), other_slot.to_variable_mat());

                // add the message into the other variable
                var other_variable = &self.variables.items[other_slot.variable_id];
                blasfeo.vecad(1.0, other_slot.to_variable_vec(), other_variable.info_vec());
                blasfeo.gead(1.0, other_slot.to_variable_mat(), other_variable.info_mat());

                // convert from info form to covariance form
                var other_mean = other_variable.mean();
                var other_cov = other_variable.cov();
                try infoToCov(other_variable.info_vec(),
                              other_variable.info_mat(),
                              other_mean,
                              other_cov);
                
                // // compute the KL divergence between previous
                // // https://mr-easy.github.io/2020-04-16-kl-divergence-between-2-gaussian-distributions/
                // // p = prev, q = current
                // var delta = try ARENA.vector(other_variable.dim);
                // var infoq_delta = try ARENA.vector(other_variable.dim);
                // var infoq_sigmap = try ARENA.matrix(other_variable.dim, other_variable.dim);

                // blasfeo.axpy(-1.0, msg_mean, other_variable.mean(), delta);
                // blasfeo.gemv_n(1.0, other_variable.info_mat(), delta, 0.0, delta, infoq_delta);
                // var deltat_infoq_delta = blasfeo.dot(delta, infoq_delta);

                // blasfeo.gemm_nn(1.0, other_variable.info_mat(), other_variable.cov(), 0.0, other_variable.cov(), infoq_sigmap);
                // const tr_infoq_sigmap = blasfeo.trace(infoq_sigmap);
                // const kl = 0.5 * (std.math.ln(msg_cov_det / variable.cov_det) - @intToFloat(f64, other_variable.dim) + deltat_infoq_delta + tr_infoq_sigmap);
            }
        }
    }
};

test "add variable" {
    var graph = FactorGraph.init(std.heap.c_allocator);
    defer graph.deinit();

    _ = try graph.addVariable(4);
    try std.testing.expect(graph.variables.items.len == 1);
}

test "add variables" {
    var graph = FactorGraph.init(std.heap.c_allocator);
    defer graph.deinit();

    _ = try graph.addVariable(4);
    _ = try graph.addVariable(8);
    _ = try graph.addVariable(4);
    try std.testing.expect(graph.variables.items.len == 3);
}

test "add factors" {
    var graph = FactorGraph.init(std.heap.c_allocator);
    defer graph.deinit();

    const v0 = try graph.addVariable(4);
    const v1 = try graph.addVariable(8);
    const v2 = try graph.addVariable(4);

    const f012 = try graph.addFactor(&.{ v0.id, v1.id, v2.id });
    const f0 = try graph.addFactor(&.{v0.id});
    const f2 = try graph.addFactor(&.{v2.id});

    _ = f012;
    _ = f0;
    _ = f2;
}

test "variable factors iteration" {
    var graph = FactorGraph.init(std.heap.c_allocator);
    defer graph.deinit();

    const v0 = try graph.addVariable(4);
    const v1 = try graph.addVariable(8);
    const v2 = try graph.addVariable(4);

    const f012 = try graph.addFactor(&.{ v0.id, v1.id, v2.id });
    const f0 = try graph.addFactor(&.{v0.id});
    const f12 = try graph.addFactor(&.{ v1.id, v2.id });

    {
        var it = graph.variableFactorsIterator(v0.id);
        var expected_factor_ids: [2]usize = .{ f012.id, f0.id };
        for (expected_factor_ids) |factor_id| {
            try std.testing.expect(it.next().?.id == factor_id);
        }
    }
    {
        var it = graph.variableFactorsIterator(v1.id);
        var expected_factor_ids: [2]usize = .{ f012.id, f12.id };
        for (expected_factor_ids) |factor_id| {
            try std.testing.expect(it.next().?.id == factor_id);
        }
    }
    {
        var it = graph.variableFactorsIterator(v2.id);
        var expected_factor_ids: [2]usize = .{ f012.id, f12.id };
        for (expected_factor_ids) |factor_id| {
            try std.testing.expect(it.next().?.id == factor_id);
        }
    }
}

test "variable init" {
    var graph = FactorGraph.init(std.heap.c_allocator);
    defer graph.deinit();

    const v0 = try graph.addVariable(2);
    // zig fmt: off
    const mean = [_]f64 { 0.1, 1 };
    const cov = [_]f64 {
        0.50, 0.25,
        0.25, 0.50
    };
    // zig fmt: on

    try v0.init(mean[0..], cov[0..]);

    // test mean and cov set properly
    var vmean = v0.mean();
    {
        var i: c_int = 0;
        while (i < 2) : (i += 1) {
            try std.testing.expect(mean[@intCast(usize, i)] == blasfeo.vecel(vmean, i));
        }
    }
    var vcov = v0.cov();
    {
        var r: c_int = 0;
        while (r < 2) : (r += 1) {
            var c: c_int = 0;
            while (c < 2) : (c += 1) {
                try std.testing.expect(cov[@intCast(usize, 2 * r + c)] == blasfeo.matel(vcov, r, c));
            }
        }
    }

    // check prior_info_mat * mean = prior_vec
    var work_vec = try ARENA.vector(v0.dim);
    blasfeo.gemv_n(1.0, v0.prior_info_mat(), v0.mean(), 0.0, v0.mean(), work_vec);
    try std.testing.expect(blasfeo.isVecEq(work_vec, v0.prior_info_vec(), 1e-10));

    // check info = prior_info
    try std.testing.expect(blasfeo.isVecEq(v0.prior_info_vec(), v0.info_vec(), 0.0));
    try std.testing.expect(blasfeo.isMatEq(v0.prior_info_mat(), v0.info_mat(), 0.0));
}

// test "propagate variable" {
//     var graph = FactorGraph.init(std.heap.c_allocator);
//     defer graph.deinit();

//     const v0 = try graph.addVariable(2);
//     // zig fmt: off
//     const mean = [_]f64 { 1.0, 1.0 };
//     const cov = [_]f64 {
//         0.50, 0.00,
//         0.00, 0.50
//     };

//     try v0.init(mean[0..], cov[0..]);

//     // add a factor that brings this variable towards zero
//     // x - r0 ~ N(0, factor_sigma)
//     // log p(x) ~ || x - r0 ||²_Q
//     //          ~ || Q½ (x - r0 ) ||²
//     //          ~ x.t Q x - 2 Q r0
//     // => info_mat = Q, info_vec = Q r0

//     var factor = try graph.addFactor(&.{v0.id});
//     const factor_info_vec = [_]f64{0.0, 0.0};
//     const factor_info_mat = [_]f64{
//         100, 0.0,
//         0.0, 100
//     };

//     factor.initWithSlice(factor_info_vec[0..], factor_info_mat[0..]);
//     // zig fmt: on

//     try graph.propagateVariable(v0.id);

//     std.debug.print("New mean\n", .{});
//     blasfeo.print_vec(v0.mean());
//     std.debug.print("New cov\n", .{});
//     blasfeo.print_mat(v0.cov());
// }

test "propagate two variables" {
    var graph = FactorGraph.init(std.heap.c_allocator);
    defer graph.deinit();

    const v0 = try graph.addVariable(2);
    const v1 = try graph.addVariable(2);
    try std.testing.expect(v0.id != v1.id);
    // zig fmt: off

    {
        const mean = [_]f64 { 1.0, 0.0 };
        const cov = [_]f64 {
            0.01,  0.0,
            0.0, 100.0
        };
        try v0.init(mean[0..], cov[0..]);
    }

    {
        const mean = [_]f64 { 0.0, 1.0 };
        const cov = [_]f64 {
            100.0, 0.00,
             0.0,  0.01
        };
        try v1.init(mean[0..], cov[0..]);
    }

    // add a factor that brings these variables equal
    // log p(x1,x2) ~ || x1 - x2 ||_Q
    //              ~ || J (x1, x2) ||_Q
    //              ~ (x1, x2) JtJ (x1, x2)
    // where J = [ I, -I ]

    var factor = try graph.addFactor(&.{v0.id, v1.id});    

    try ARENA.pushEnv();
    defer ARENA.popEnv();

    var J = try ARENA.matrix(2, 4);
    blasfeo.gese(0.0, J);
    blasfeo.diare(1.0, J.block(0, 0, 2, 2));
    blasfeo.diare(-1.0, J.block(0, 2, 2, 2));
    std.debug.print("J\n", .{});
    blasfeo.print_mat(J);
    
    var info_mat = try ARENA.matrix(4, 4);
    var zero44 = try ARENA.matrix(4,4);
    blasfeo.gese(0.0, zero44);
    blasfeo.gemm_tn(1.0, J, J, 0.0, zero44, info_mat);
    var info_vec = try ARENA.vector(4);
    blasfeo.vecse(0.0, info_vec);

    std.debug.print("Info mat\n", .{});
    blasfeo.print_mat(info_mat);
    std.debug.print("Info vec\n", .{});
    blasfeo.print_vec(info_vec);
    
    factor.relinearize(info_vec, info_mat);
    std.debug.print("Factor info initialized to:\n", .{});
    blasfeo.print_mat(factor.info_mat());
    // zig fmt: on

    var it: usize = 0;
    while (it < 1) : (it += 1) {
        try graph.propagateVariable(v0.id);
        try graph.propagateVariable(v1.id);
        std.debug.print("New v0 mean\n", .{});
        blasfeo.print_vec(v0.mean());
        std.debug.print("New v1 mean\n", .{});
        blasfeo.print_vec(v1.mean());
    }

    // std.debug.print("New v0 cov\n", .{});
    // blasfeo.print_mat(v0.cov());
}

// test "allocVecs" {
//     var dvec1: blasfeo.Dvec = undefined;
//     var dvec2: blasfeo.Dvec = undefined;
//     var dvecs: [2]*blasfeo.Dvec = .{ &dvec1, &dvec2 };
//     try allocVecs(std.heap.c_allocator, 4, &dvecs);
//     for (dvecs) |dvec| {
//         try std.testing.expect(dvec.mem != 0);
//         if (dvec.mem != 0) {
//             std.heap.c_allocator.destroy(dvec.mem);
//         }
//     }
// }
