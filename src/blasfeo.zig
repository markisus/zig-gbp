const std = @import("std");
const cblasfeo = @cImport({
    // todo: coordinate these values
    @cInclude("blasfeo.h");
    @cInclude("blasfeo_d_blasfeo_ref_api.h");
});

pub const CACHE_LINE_SIZE = cblasfeo.CACHE_LINE_SIZE;

pub const Dmat = cblasfeo.blasfeo_dmat;
pub const Dvec = cblasfeo.blasfeo_dvec;

pub fn removeConst(comptime T: type, ptr: *const T) *T {
    return @intToPtr(*T, @ptrToInt(ptr));
}

// A view of a blasfeo matrix
pub const Matrix = struct {
    const Self = @This();

    dmat: *cblasfeo.blasfeo_dmat = undefined,
    row_start: c_int = 0,
    col_start: c_int = 0,
    rows: c_int = 0,
    cols: c_int = 0,

    pub fn attach(dmat: *Dmat) Matrix {
        return Matrix{ .dmat = dmat, .rows = dmat.m, .cols = dmat.n };
    }

    pub fn block(self: *Self, row_start_1: c_int, col_start_1: c_int, rows_1: c_int, cols_1: c_int) Matrix {
        return .{ .dmat = self.dmat, .row_start = self.row_start + row_start_1, .col_start = self.col_start + col_start_1, .rows = rows_1, .cols = cols_1 };
    }

    pub fn setFromSlice(self: *Self, data: []const f64) void {
        std.debug.assert(data.len == self.rows * self.cols);
        pack_tran_dmat(removeConst(f64, &data.ptr[0]), self.cols, self.*);
    }
};

// A view of blasfeo vector
pub const Vector = struct {
    const Self = @This();

    dvec: *cblasfeo.blasfeo_dvec,
    len: c_int = 0,
    start: c_int = 0,

    pub fn attach(dvec: *Dvec) Vector {
        return Vector{ .dvec = dvec, .len = dvec.m };
    }

    pub fn segment(self: *const Self, start_1: c_int, len_1: c_int) Vector {
        std.debug.assert(0 <= start_1 and start_1 + len_1 <= self.len);
        return .{ .dvec = self.dvec, .len = len_1, .start = self.start + start_1 };
    }
    pub fn head(self: *const Self, len_1: c_int) Vector {
        std.debug.assert(len_1 <= self.len);
        return .{ .dvec = self.dvec, .len = len_1, .start = self.start };
    }
    pub fn tail(self: *const Self, len_1: c_int) Vector {
        std.debug.assert(len_1 <= self.len);
        return .{ .dvec = self.dvec, .len = len_1, .start = self.start + self.len - len_1 };
    }

    pub fn setFromSlice(self: *Self, data: []const f64) void {
        std.debug.assert(data.len == self.len);
        pack_dvec(removeConst(f64, &data.ptr[0]), 1, self.*);
    }
};

pub fn matel(matrix: Matrix, ai_: c_int, aj_: c_int) f64 {
    // manual translation of DMATEL under MF_PANELMAJ & LA_HIGH_PERFORMANCE
    // #define BLASFEO_DMATEL(sA,ai,aj) ((sA)->pA[((ai)-((ai)&(D_PS-1)))*(sA)->cn+(aj)*D_PS+((ai)&(D_PS-1))])
    // this macro does not get translated correctly because index must be usize
    var sA = matrix.dmat;
    const D_PS = cblasfeo.D_PS;
    const ai = ai_ + matrix.row_start;
    const aj = aj_ + matrix.col_start;
    return (sA.pA[@intCast(usize, ((ai) - ((ai) & (D_PS - 1))) * sA.cn + (aj) * D_PS + ((ai) & (D_PS - 1)))]);
}

pub fn print_mat(matrix: Matrix) void {
    cblasfeo.blasfeo_print_dmat(matrix.rows, matrix.cols, matrix.dmat, matrix.row_start, matrix.col_start);
}

pub fn print_vec(vector: Vector) void {
    cblasfeo.blasfeo_print_dvec(vector.len, vector.dvec, vector.start);
}

pub fn vecel(vector: Vector, ai_: c_int) f64 {
    // manual translation of DMATEL under MF_PANELMAJ & LA_HIGH_PERFORMANCE
    // #define BLASFEO_DVECEL(sa,ai) ((sa)->pa[ai])
    // this macro does not get translated correctly because index must be usize
    var sa = vector.dvec;
    const ai = ai_ + vector.start;
    return sa.pa[@intCast(usize, ai)];
}

pub fn setId(m: Matrix) void {
    gese(0.0, m);
    diare(1.0, m);
}

pub fn isMatEq(m1: Matrix, m2: Matrix, tol: f64) bool {
    if (m1.rows != m2.rows) return false;
    if (m1.cols != m2.cols) return false;

    var c: c_int = 0;
    while (c < m1.cols) : (c += 1) {
        var r: c_int = 0;
        while (r < m1.rows) : (r += 1) {
            if (@fabs(matel(m1, r, c) - matel(m2, r, c)) > tol) return false;
        }
    }
    return true;
}

pub fn isVecEq(v1: Vector, v2: Vector, tol: f64) bool {
    if (v1.len != v2.len) return false;

    var c: c_int = 0;
    while (c < v1.len) : (c += 1) {
        if (@fabs(vecel(v1, c) - vecel(v2, c)) > tol) return false;
    }
    return true;
}

pub fn trace(matrix: Matrix) f64 {
    var result: f64 = 0;
    var idx: c_int = 0;
    while (idx < @minimum(matrix.rows, matrix.cols)) : (idx += 1) {
        result += matel(matrix, idx, idx);
    }
    return result;
}

pub fn diagonalProduct(matrix: Matrix) f64 {
    var result: f64 = 1;
    var idx: c_int = 0;
    while (idx < @minimum(matrix.rows, matrix.cols)) : (idx += 1) {
        result *= matel(matrix, idx, idx);
    }
    return result;
}

pub fn memsize_mat(m: c_int, n: c_int) usize {
    return cblasfeo.blasfeo_memsize_dmat(m, n);
}
// returns the memory size (in bytes) needed for a dvec
pub fn memsize_vec(m: c_int) usize {
    return cblasfeo.blasfeo_memsize_dvec(m);
}

pub fn pack_tran_dmat(A: *f64, lda: c_int, B: Matrix) void {
    // todo: fixme, is something wrong about order of rows and cols?
    cblasfeo.blasfeo_pack_tran_dmat(B.rows, B.cols, A, lda, B.dmat, B.row_start, B.col_start);
}

pub fn pack_dvec(x: *f64, xi: c_int, y: Vector) void {
    cblasfeo.blasfeo_pack_dvec(y.len, x, xi, y.dvec, y.start);
}
pub fn gese(alpha: f64, A: Matrix) void {
    cblasfeo.blasfeo_dgese(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start);
}
pub fn diare(alpha: f64, A: Matrix) void {
    cblasfeo.blasfeo_ddiare(A.rows, alpha, A.dmat, A.row_start, A.col_start);
}
pub fn create_dmat(m: c_int, n: c_int, dmat: *Dmat, memory: *anyopaque) void {
    cblasfeo.blasfeo_create_dmat(m, n, dmat, memory);
}
pub fn create_dvec(m: c_int, dvec: *Dvec, memory: *anyopaque) void {
    cblasfeo.blasfeo_create_dvec(m, dvec, memory);
}
pub fn vecse(alpha: f64, x: Vector) void {
    cblasfeo.blasfeo_dvecse(x.len, alpha, x.dvec, x.start);
}
pub fn veccp(x: Vector, y: Vector) void {
    cblasfeo.blasfeo_dveccp(x.len, x.dvec, x.start, y.dvec, y.start);
}
pub fn gecp(A: Matrix, B: Matrix) void {
    cblasfeo.blasfeo_dgecp(A.rows, A.cols, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
}
pub fn potrf_l(C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dpotrf_l(C.rows, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_llnn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_llnn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_lltn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_lltn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsv_lnn(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrsv_lnn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trsv_ltn(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrsv_ltn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn gemv_n(alpha: f64, A: Matrix, x: Vector, beta: f64, y: Vector, z: Vector) void {
    cblasfeo.blasfeo_dgemv_n(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
}
pub fn gemm_tn(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
    // todo: Is this right?
    cblasfeo.blasfeo_dgemm_tn(D.rows, D.cols, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn gead(alpha: f64, A: Matrix, B: Matrix) void {
    cblasfeo.blasfeo_dgead(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
}
pub fn vecad(alpha: f64, x: Vector, y: Vector) void {
    cblasfeo.blasfeo_dvecad(x.len, alpha, x.dvec, x.start, y.dvec, y.start);
}
pub fn axpy(alpha: f64, x: Vector, y: Vector, z: Vector) void {
    cblasfeo.blasfeo_daxpy(x.len, alpha, x.dvec, x.start, y.dvec, y.start, z.dvec, z.start);
}
pub fn dot(x: Vector, y: Vector) f64 {
    return cblasfeo.blasfeo_ddot(x.len, x.dvec, x.start, y.dvec, y.start);
}
pub fn gemm_nn(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dgemm_nn(A.rows, A.cols, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}

// pub fn axpby(alpha: f64, x: Vector, beta: f64, y: Vector, z: Vector) void {
//     cblasfeo.blasfeo_daxpby(x.len, alpha, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
// }
// pub fn vecmul(x: Vector, y: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dvecmul(x.len, x.dvec, x.start, y.dvec, y.start, z.dvec, z.start);
// }
// pub fn vecmulacc(x: Vector, y: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dvecmulacc(x.len, x.dvec, x.start, y.dvec, y.start, z.dvec, z.start);
// }
// pub fn vecmuldot(x: Vector, y: Vector, z: Vector) f64 {
//     cblasfeo.blasfeo_dvecmuldot(x.len, x.dvec, x.start, y.dvec, y.start, z.dvec, z.start);
// }
// pub fn rotg(a: f64, b: f64, c: *f64, s: *f64) void {
//     cblasfeo.blasfeo_drotg(a, b, c, s);
// }
// pub fn colrot(A: Matrix, aj0: c_int, aj1: c_int, c: f64, s: f64) void {
//     cblasfeo.blasfeo_dcolrot(A.rows, A.dmat, A.row_start, aj0, aj1, c, s);
// }
// pub fn rowrot(A: Matrix, ai0: c_int, ai1: c_int, c: f64, s: f64) void {
//     cblasfeo.blasfeo_drowrot(A.rows, A.dmat, ai0, ai1, A.col_start, c, s);
// }
// pub fn gemv_t(alpha: f64, A: Matrix, x: Vector, beta: f64, y: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dgemv_t(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
// }
// pub fn trsv_lnn_mn(A: Matrix, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dtrsv_lnn_mn(A.rows, A.cols, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn trsv_ltn_mn(A: Matrix, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dtrsv_ltn_mn(A.rows, A.cols, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn trsv_lnu(A: Matrix, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dtrsv_lnu(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn trsv_ltu(A: Matrix, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dtrsv_ltu(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn trsv_unn(A: Matrix, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dtrsv_unn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn trsv_utn(A: Matrix, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dtrsv_utn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn trmv_lnn(A: Matrix, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dtrmv_lnn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn trmv_lnu(A: Matrix, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dtrmv_lnu(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn trmv_ltn(A: Matrix, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dtrmv_ltn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn trmv_ltu(A: Matrix, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dtrmv_ltu(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn trmv_unn(A: Matrix, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dtrmv_unn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn trmv_utn(A: Matrix, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dtrmv_utn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn gemv_nt(alpha_n: f64, alpha_t: f64, A: Matrix, x_n: Vector, x_t: Vector, beta_n: f64, beta_t: f64, y_n: Vector, y_t: Vector, z_n: Vector, z_t: Vector) void {
//     cblasfeo.blasfeo_dgemv_nt(A.rows, A.cols, alpha_n, alpha_t, A.dmat, A.row_start, A.col_start, x_n.dvec, x_n.start, x_t.dvec, x_t.start, beta_n, beta_t, y_n.dvec, y_n.start, y_t.dvec, y_t.start, z_n.dvec, z_n.start, z_t.dvec, z_t.start);
// }
// pub fn symv_l(alpha: f64, A: Matrix, x: Vector, beta: f64, y: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dsymv_l(A.rows, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
// }
// pub fn symv_l_mn(alpha: f64, A: Matrix, x: Vector, beta: f64, y: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dsymv_l_mn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
// }
// pub fn symv_u(alpha: f64, A: Matrix, x: Vector, beta: f64, y: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dsymv_u(A.rows, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
// }
// pub fn ger(alpha: f64, x: Vector, y: Vector, C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dger(C.rows, C.cols, alpha, x.dvec, x.start, y.dvec, y.start, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn gemv_d(alpha: f64, A: Vector, ai: c_int, x: Vector, beta: f64, y: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dgemv_d(A.len, alpha, A.dvec, ai, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
// }
// pub fn gemm_nt(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dgemm_nt(A.rows, A.cols, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn gemm_tt(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dgemm_tt(A.rows, A.cols, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn syrk_ln(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dsyrk_ln(A.rows, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn syrk_ln_mn(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dsyrk_ln_mn(A.rows, A.cols, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn syrk_lt(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dsyrk_lt(A.rows, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn syrk_un(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dsyrk_un(A.rows, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn syrk_ut(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dsyrk_ut(A.rows, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_llnn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_llnn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_llnu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_llnu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_lltn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_lltn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_lltu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_lltu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_lunn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_lunn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_lunu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_lunu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_lutn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_lutn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_lutu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_lutu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_rlnn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_rlnn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_rlnu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_rlnu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_rltn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_rltn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_rltu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_rltu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_runn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_runn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_runu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_runu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_rutn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_rutn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trmm_rutu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrmm_rutu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_llnu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_llnu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_lltu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_lltu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_lunn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_lunn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_lunu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_lunu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_lutn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_lutn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_lutu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_lutu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_rlnn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_rlnn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_rlnu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_rlnu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_rltn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_rltn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_rltu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_rltu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_runn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_runn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_runu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_runu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_rutn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_rutn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn trsm_rutu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dtrsm_rutu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn gemm_dn(alpha: f64, A: Vector, ai: c_int, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dgemm_dn(B.rows, B.cols, alpha, A.dvec, ai, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn gemm_nd(alpha: f64, A: Matrix, B: Vector, bi: c_int, beta: f64, C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dgemm_nd(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dvec, bi, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn potrf_l_mn(C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dpotrf_l_mn(C.rows, C.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn potrf_u(C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dpotrf_u(C.rows, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn syrk_dpotrf_ln(A: Matrix, B: Matrix, C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dsyrk_dpotrf_ln(A.rows, D.cols, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn syrk_dpotrf_ln_mn(A: Matrix, B: Matrix, C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dsyrk_dpotrf_ln_mn(A.rows, A.cols, D.cols, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn getrf_np(C: Matrix, D: Matrix) void {
//     cblasfeo.blasfeo_dgetrf_np(C.rows, C.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
// }
// pub fn getrf_rp(C: Matrix, D: Matrix, ipiv: *c_int) void {
//     cblasfeo.blasfeo_dgetrf_rp(C.rows, C.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start, ipiv);
// }
// pub fn geqrf(C: Matrix, D: Matrix, work: *anyopaque) void {
//     cblasfeo.blasfeo_dgeqrf(C.rows, C.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start, work);
// }
// pub fn orglq(C: Matrix, D: Matrix, work: *anyopaque) void {
//     cblasfeo.blasfeo_dorglq(C.rows, C.cols, D.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start, work);
// }
// pub fn gelqf(C: Matrix, D: Matrix, work: *anyopaque) void {
//     cblasfeo.blasfeo_dgelqf(C.rows, C.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start, work);
// }
// pub fn gelqf_pd(C: Matrix, D: Matrix, work: *anyopaque) void {
//     cblasfeo.blasfeo_dgelqf_pd(C.rows, C.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start, work);
// }
// pub fn gelqf_pd_la(n1: c_int, L: Matrix, A: Matrix, work: *anyopaque) void {
//     cblasfeo.blasfeo_dgelqf_pd_la(L.rows, n1, L.dmat, L.row_start, L.col_start, A.dmat, A.row_start, A.col_start, work);
// }
// pub fn gelqf_pd_lla(n1: c_int, L0: Matrix, L1: Matrix, A: Matrix, work: *anyopaque) void {
//     cblasfeo.blasfeo_dgelqf_pd_lla(L0.rows, n1, L0.dmat, L0.row_start, L0.col_start, L1.dmat, L1.row_start, L1.col_start, A.dmat, A.row_start, A.col_start, work);
// }

// pub fn pack_dmat(A: *f64, lda: c_int, B: Matrix) void {
//     cblasfeo.blasfeo_pack_dmat(B.rows, B.cols, A, lda, B.dmat, B.row_start, B.col_start);
// }
// pub fn pack_l_dmat(A: *f64, lda: c_int, B: Matrix) void {
//     cblasfeo.blasfeo_pack_l_dmat(B.rows, B.cols, A, lda, B.dmat, B.row_start, B.col_start);
// }
// pub fn pack_u_dmat(A: *f64, lda: c_int, B: Matrix) void {
//     cblasfeo.blasfeo_pack_u_dmat(B.rows, B.cols, A, lda, B.dmat, B.row_start, B.col_start);
// }
// pub fn unpack_dmat(A: Matrix, B: *f64, ldb: c_int) void {
//     cblasfeo.blasfeo_unpack_dmat(A.rows, A.cols, A.dmat, A.row_start, A.col_start, B, ldb);
// }
// pub fn unpack_tran_dmat(A: Matrix, B: *f64, ldb: c_int) void {
//     cblasfeo.blasfeo_unpack_tran_dmat(A.rows, A.cols, A.dmat, A.row_start, A.col_start, B, ldb);
// }
// pub fn unpack_dvec(x: Vector, y: *f64, yi: c_int) void {
//     cblasfeo.blasfeo_unpack_dvec(x.len, x.dvec, x.start, y, yi);
// }
// pub fn gesc(alpha: f64, A: Matrix) void {
//     cblasfeo.blasfeo_dgesc(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start);
// }
// pub fn gecpsc(alpha: f64, A: Matrix, B: Matrix) void {
//     cblasfeo.blasfeo_dgecpsc(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
// }
// pub fn trcp_l(A: Matrix, B: Matrix) void {
//     cblasfeo.blasfeo_dtrcp_l(A.rows, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
// }
// pub fn trcpsc_l(alpha: f64, A: Matrix, B: Matrix) void {
//     cblasfeo.blasfeo_dtrcpsc_l(A.rows, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
// }
// pub fn trsc_l(alpha: f64, A: Matrix) void {
//     cblasfeo.blasfeo_dtrsc_l(A.rows, alpha, A.dmat, A.row_start, A.col_start);
// }
// pub fn getr(A: Matrix, B: Matrix) void {
//     cblasfeo.blasfeo_dgetr(A.rows, A.cols, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
// }
// pub fn trtr_l(A: Matrix, B: Matrix) void {
//     cblasfeo.blasfeo_dtrtr_l(A.rows, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
// }
// pub fn trtr_u(A: Matrix, B: Matrix) void {
//     cblasfeo.blasfeo_dtrtr_u(A.rows, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
// }
// pub fn diain(alpha: f64, x: Vector, A: Matrix) void {
//     cblasfeo.blasfeo_ddiain(A.rows, alpha, x.dvec, x.start, A.dmat, A.row_start, A.col_start);
// }
// pub fn diain_sp(alpha: f64, x: Vector, idx: *c_int, D: Matrix) void {
//     cblasfeo.blasfeo_ddiain_sp(D.rows, alpha, x.dvec, x.start, idx, D.dmat, D.row_start, D.col_start);
// }
// pub fn diaex(alpha: f64, A: Matrix, x: Vector) void {
//     cblasfeo.blasfeo_ddiaex(A.rows, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start);
// }
// pub fn diaex_sp(alpha: f64, idx: *c_int, D: Matrix, x: Vector) void {
//     cblasfeo.blasfeo_ddiaex_sp(D.rows, alpha, idx, D.dmat, D.row_start, D.col_start, x.dvec, x.start);
// }
// pub fn diaad(alpha: f64, x: Vector, A: Matrix) void {
//     cblasfeo.blasfeo_ddiaad(A.rows, alpha, x.dvec, x.start, A.dmat, A.row_start, A.col_start);
// }
// pub fn diaad_sp(alpha: f64, x: Vector, idx: *c_int, D: Matrix) void {
//     cblasfeo.blasfeo_ddiaad_sp(D.rows, alpha, x.dvec, x.start, idx, D.dmat, D.row_start, D.col_start);
// }
// pub fn diaadin_sp(alpha: f64, x: Vector, y: Vector, idx: *c_int, D: Matrix) void {
//     cblasfeo.blasfeo_ddiaadin_sp(D.rows, alpha, x.dvec, x.start, y.dvec, y.start, idx, D.dmat, D.row_start, D.col_start);
// }
// pub fn rowin(alpha: f64, x: Vector, A: Matrix) void {
//     cblasfeo.blasfeo_drowin(A.rows, alpha, x.dvec, x.start, A.dmat, A.row_start, A.col_start);
// }
// pub fn rowex(alpha: f64, A: Matrix, x: Vector) void {
//     cblasfeo.blasfeo_drowex(A.rows, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start);
// }
// pub fn rowad(alpha: f64, x: Vector, A: Matrix) void {
//     cblasfeo.blasfeo_drowad(A.rows, alpha, x.dvec, x.start, A.dmat, A.row_start, A.col_start);
// }
// pub fn rowad_sp(alpha: f64, x: Vector, idx: *c_int, D: Matrix) void {
//     cblasfeo.blasfeo_drowad_sp(D.rows, alpha, x.dvec, x.start, idx, D.dmat, D.row_start, D.col_start);
// }
// pub fn rowsw(A: Matrix, C: Matrix) void {
//     cblasfeo.blasfeo_drowsw(A.rows, A.dmat, A.row_start, A.col_start, C.dmat, C.row_start, C.col_start);
// }
// pub fn rowpe(ipiv: *c_int, A: Matrix) void {
//     cblasfeo.blasfeo_drowpe(A.rows, ipiv, A.dmat);
// }
// pub fn rowpei(ipiv: *c_int, A: Matrix) void {
//     cblasfeo.blasfeo_drowpei(A.rows, ipiv, A.dmat);
// }
// pub fn colex(A: Matrix, x: Vector) void {
//     cblasfeo.blasfeo_dcolex(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start);
// }
// pub fn colin(x: Vector, A: Matrix) void {
//     cblasfeo.blasfeo_dcolin(A.rows, x.dvec, x.start, A.dmat, A.row_start, A.col_start);
// }
// pub fn colad(alpha: f64, x: Vector, A: Matrix) void {
//     cblasfeo.blasfeo_dcolad(A.rows, alpha, x.dvec, x.start, A.dmat, A.row_start, A.col_start);
// }
// pub fn colsc(alpha: f64, A: Matrix) void {
//     cblasfeo.blasfeo_dcolsc(A.rows, alpha, A.dmat, A.row_start, A.col_start);
// }
// pub fn colsw(A: Matrix, C: Matrix) void {
//     cblasfeo.blasfeo_dcolsw(A.rows, A.dmat, A.row_start, A.col_start, C.dmat, C.row_start, C.col_start);
// }
// pub fn colpe(ipiv: *c_int, A: Matrix) void {
//     cblasfeo.blasfeo_dcolpe(A.rows, ipiv, A.dmat);
// }
// pub fn colpei(ipiv: *c_int, A: Matrix) void {
//     cblasfeo.blasfeo_dcolpei(A.rows, ipiv, A.dmat);
// }
// pub fn vecsc(alpha: f64, x: Vector) void {
//     cblasfeo.blasfeo_dvecsc(x.len, alpha, x.dvec, x.start);
// }
// pub fn veccpsc(alpha: f64, x: Vector, y: Vector) void {
//     cblasfeo.blasfeo_dveccpsc(x.len, alpha, x.dvec, x.start, y.dvec, y.start);
// }
// pub fn vecad_sp(alpha: f64, x: Vector, idx: *c_int, z: Vector) void {
//     cblasfeo.blasfeo_dvecad_sp(x.len, alpha, x.dvec, x.start, idx, z.dvec, z.start);
// }
// pub fn vecin_sp(alpha: f64, x: Vector, idx: *c_int, z: Vector) void {
//     cblasfeo.blasfeo_dvecin_sp(x.len, alpha, x.dvec, x.start, idx, z.dvec, z.start);
// }
// pub fn vecex_sp(alpha: f64, idx: *c_int, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dvecex_sp(x.len, alpha, idx, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn vecexad_sp(alpha: f64, idx: *c_int, x: Vector, z: Vector) void {
//     cblasfeo.blasfeo_dvecexad_sp(x.len, alpha, idx, x.dvec, x.start, z.dvec, z.start);
// }
// pub fn vecze(m: Vector, v: Vector, e: Vector) void {
//     cblasfeo.blasfeo_dvecze(m.len, m.dvec, m.start, v.dvec, v.start, e.dvec, e.start);
// }
// pub fn vecnrm_inf(x: Vector, ptr_norm: *f64) void {
//     cblasfeo.blasfeo_dvecnrm_inf(x.len, x.dvec, x.start, ptr_norm);
// }
// pub fn vecpe(ipiv: *c_int, x: Vector) void {
//     cblasfeo.blasfeo_dvecpe(x.len, ipiv, x.dvec, x.start);
// }
// pub fn vecpei(ipiv: *c_int, x: Vector) void {
//     cblasfeo.blasfeo_dvecpei(x.len, ipiv, x.dvec, x.start);
// }

pub const Arena = struct {
    const Self = @This();
    const UNINTIALIZED: [0]u8 = .{};
    const ARENA_SIZE_KB = 64;

    env_stack: std.ArrayList(usize),
    allocator: std.mem.Allocator,
    bytes: []u8 align(CACHE_LINE_SIZE) = UNINTIALIZED[0..],
    bytes_used: usize = 0,
    bytes_allocated: bool = false,

    pub fn init(allocator: std.mem.Allocator) Self {
        var result = Self{ .env_stack = std.ArrayList(usize).init(allocator), .allocator = allocator };
        return result;
    }

    // bump allocate a number of bytes, aligned to cache
    pub fn alloc(self: *Self, num_bytes: usize) ![]align(CACHE_LINE_SIZE) u8 {
        if (!self.bytes_allocated) {
            self.bytes = try self.allocator.alignedAlloc(u8, CACHE_LINE_SIZE, 1024 * ARENA_SIZE_KB);
            self.bytes_allocated = true;
        }

        // round up to multiple of cache line size
        var num_bytes_aligned = @divFloor(num_bytes, CACHE_LINE_SIZE) * CACHE_LINE_SIZE;

        if (num_bytes % CACHE_LINE_SIZE != 0) {
            num_bytes_aligned += CACHE_LINE_SIZE;
        }

        std.debug.assert(num_bytes_aligned != 0);

        const start_byte = self.bytes_used;
        self.bytes_used += num_bytes_aligned;
        if (self.bytes_used >= self.bytes.len) {
            return error.OutOfMemory;
        }

        // assert alignment to cache line
        std.debug.assert(@ptrToInt(&self.bytes[start_byte]) % CACHE_LINE_SIZE == 0);

        var result: []align(CACHE_LINE_SIZE) u8 = @alignCast(CACHE_LINE_SIZE, self.bytes[start_byte..]);
        result.len = num_bytes;

        return result;
    }

    pub fn pushEnv(self: *Self) !void {
        try self.env_stack.append(self.bytes_used);
    }

    pub fn popEnv(self: *Self) void {
        std.debug.assert(self.env_stack.items.len >= 1);
        const last_env = self.env_stack.pop();
        self.bytes_used = last_env;
    }

    pub fn deinit(self: *Self) void {
        self.env_stack.deinit();
        if (self.bytes_allocated) {
            self.allocator.free(self.bytes);
        }
    }

    pub fn matrix(self: *Self, rows: c_int, cols: c_int) !Matrix {
        const checkpoint = self.bytes_used;
        errdefer self.bytes_used = checkpoint;

        const header_memsize: usize = @sizeOf(Dmat);
        const storage_memsize: usize = memsize_mat(rows, cols);
        const header_bytes = try self.alloc(header_memsize);
        const storage_bytes = try self.alloc(storage_memsize);
        create_dmat(rows, cols, @ptrCast(*Dmat, header_bytes.ptr), storage_bytes.ptr);
        return Matrix.attach(@ptrCast(*Dmat, header_bytes.ptr));
    }

    pub fn vector(self: *Self, len: c_int) !Vector {
        const checkpoint = self.bytes_used;
        errdefer self.bytes_used = checkpoint;

        const header_memsize: usize = @sizeOf(Dvec);
        const storage_memsize: usize = memsize_vec(len);
        const header_bytes = try self.alloc(header_memsize);
        const storage_bytes = try self.alloc(storage_memsize);
        create_dvec(len, @ptrCast(*Dvec, header_bytes.ptr), storage_bytes.ptr);
        return Vector.attach(@ptrCast(*Dvec, header_bytes.ptr));
    }
};
