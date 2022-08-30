const std = @import("std");
const blasfeo = @import("blasfeo.zig");

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

const FactorGraph = struct {
    const Self = @This();

    // zig fmt: off
    const Variable = struct {
        id: usize,
        mean: blasfeo.Dvec = undefined,
        cov: blasfeo.Dmat = undefined,
        info_vec: blasfeo.Dvec = undefined,
        info_mat: blasfeo.Dmat = undefined,
        prior_info_vec: blasfeo.Dvec = undefined,
        prior_info_mat: blasfeo.Dmat = undefined,
        first_factor_slot : ?usize = null,
        last_factor_slot : ?usize = null,
    };

    const Factor = struct {
        id: usize,
        info_vec : blasfeo.Dvec = undefined,
        info_mat : blasfeo.Dmat = undefined,
        slot_begin : usize = 0,
        slot_end : usize = 0
    };

    const FactorSlot = struct {
        factor_id : usize,
        variable_id : usize,
        next_same_variable_slot : ?usize = null,
        to_variable_vec : blasfeo.Dvec = undefined,
        to_variable_mat : blasfeo.Dmat = undefined,
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

    pub fn addVariable(self: *Self, dim: c_int) !usize {
        var variable = Variable{ .id = self.variables.items.len };
        var vecs: [3]*blasfeo.Dvec = .{ &variable.mean, &variable.info_vec, &variable.prior_info_vec };
        try allocVecs(self.allocator, dim, &vecs);
        errdefer {
            for (vecs) |vec| {
                self.allocator.destroy(vec.mem);
            }
        }
        var mats: [3]*blasfeo.Dmat = .{ &variable.cov, &variable.info_mat, &variable.prior_info_mat };
        try allocMats(self.allocator, dim, &mats);
        errdefer {
            for (mats) |mat| {
                self.allocator.destroy(mat.mem);
            }
        }
        try self.variables.append(variable);
        return variable.id;
    }

    pub fn addFactor(self: *Self, variable_ids: []usize) !usize {
        // std.debug.print("Adding factor with ids {d}\n", .{variable_ids});
        // factor dim is the sum of the variable dims
        var dim: c_int = 0;
        for (variable_ids) |variable_id| {
            dim += self.variables.items[variable_id].mean.m;
        }
        var factor = Factor{ .id = self.factors.items.len };

        try allocVecs(self.allocator, dim, &.{&factor.info_vec});
        errdefer {
            self.allocator.destroy(factor.info_vec.mem);
        }

        try allocMats(self.allocator, dim, &.{&factor.info_mat});
        errdefer {
            self.allocator.destroy(factor.info_mat.mem);
        }

        const slot_begin = self.factor_slots.items.len;
        errdefer {
            while (self.factor_slots.items.len > slot_begin) {
                var slot = self.factor_slots.pop();
                self.allocator.destroy(slot.to_variable_vec.mem);
                self.allocator.destroy(slot.to_variable_mat.mem);
            }
        }
        for (variable_ids) |variable_id| {
            var slot = FactorSlot{
                .factor_id = factor.id,
                .variable_id = variable_id,
            };
            const variable = &self.variables.items[variable_id];
            try allocVecs(self.allocator, variable.mean.m, &.{&slot.to_variable_vec});
            errdefer self.allocator.destroy(slot.to_variable_vec.mem);
            try allocMats(self.allocator, variable.mean.m, &.{&slot.to_variable_mat});
            errdefer self.allocator.destroy(slot.to_variable_mat.mem);
            try self.factor_slots.append(slot);
        }
        const slot_end = self.factor_slots.items.len;

        factor.slot_begin = slot_begin;
        factor.slot_end = slot_end;

        try self.factors.append(factor);

        // update internal linked lists inside the factor_slots
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

        return factor.id;
    }

    pub fn propagateVariable(self: *Self, variable_id: usize) void {
        const variable = &self.variables[variable_id];
        if (variable.first_factor_slot == null) {
            // this variable is not involved with any factors
            return;
        }

        var current_slot_idx = variable.first_factor_slot.?;
        while (current_slot_idx != variable.last_factor_slot.?) {
            var slot = &self.factor_slots.items[current_slot_idx];
            var factor = &self.factors.items[slot.factor_id];
            // compute to factor message by taking the variable info
            // and subtracting out the slot info

            // prepare a buffer for factor info

            var other_slot_idx = factor.slot_begin;
            while (other_slot_idx < factor.slot_end) : (other_slot_idx += 1) {
                // skip if other_slot_idx == current_slot_idx
                // add in the info from the slot into the factor info
            }

            // cholesky solve the factor info

            other_slot_idx = factor.slot_begin;
            while (other_slot_idx < factor.slot_end) : (other_slot_idx += 1) {
                // skip if other_slot_idx == current_slot_idx
                // subtract previous message from variable
                // add in current message to variable
                // overwite current message into the slot
                // compute kldivergence from previous, update the priority queue
            }

            if (slot.next_same_variable_slot) |next_slot_idx| {
                current_slot_idx = next_slot_idx;
            } else {
                break;
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

    const f012 = try graph.addFactor(&.{ v0, v1, v2 });
    const f0 = try graph.addFactor(&.{v0});
    const f2 = try graph.addFactor(&.{v2});

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

    const f012 = try graph.addFactor(&.{ v0, v1, v2 });
    const f0 = try graph.addFactor(&.{v0});
    const f12 = try graph.addFactor(&.{ v1, v2 });

    {
        var it = graph.variableFactorsIterator(v0);
        var expected_factor_ids: [2]usize = .{ f012, f0 };
        for (expected_factor_ids) |factor_id| {
            try std.testing.expect(it.next().?.id == factor_id);
        }
    }
    {
        var it = graph.variableFactorsIterator(v1);
        var expected_factor_ids: [2]usize = .{ f012, f12 };
        for (expected_factor_ids) |factor_id| {
            try std.testing.expect(it.next().?.id == factor_id);
        }
    }
    {
        var it = graph.variableFactorsIterator(v2);
        var expected_factor_ids: [2]usize = .{ f012, f12 };
        for (expected_factor_ids) |factor_id| {
            try std.testing.expect(it.next().?.id == factor_id);
        }
    }
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
