const std = @import("std");
const cblasfeo = @cImport({
    // todo: coordinate these values
    @cDefine("LA_HIGH_PERFORMANCE", "");
    @cDefine("MF_PANELMAJ", "");
    @cDefine("CACHE_LINE_SIZE", "64");
    @cInclude("blasfeo.h");
    @cInclude("blasfeo_d_blasfeo_ref_api.h");
});

const CACHE_LINE_SIZE = cblasfeo.CACHE_LINE_SIZE;

// A view of a blasfeo matrix
const Matrix = struct {
    const Self = @This();

    dmat: *cblasfeo.blasfeo_dmat = undefined,
    row_start: c_int = 0,
    col_start: c_int = 0,
    rows: c_int = 0,
    cols: c_int = 0,

    pub fn block(self: *Self, row_start_1: c_int, col_start_1: c_int, rows_1: c_int, cols_1: c_int) Matrix {
        return .{ .dmat = self.dmat, .row_start = self.row_start + row_start_1, .col_start = self.col_start + col_start_1, .rows = rows_1, .cols = cols_1 };
    }
};

// A view of blasfeo vector
const Vector = struct {
    const Self = @This();

    dvec: *cblasfeo.blasfeo_dvec,
    len: c_int = 0,
    start: c_int = 0,

    pub fn segment(self: *const Self, len_1: c_int, start_1: c_int) Vector {
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
};

pub fn axpy(alpha: f64, x: Vector, y: Vector, z: Vector) void {
    cblasfeo.blasfeo_daxpy(x.len, alpha, x.dvec, x.start, y.dvec, y.start, z.dvec, z.start);
}
pub fn axpby(alpha: f64, x: Vector, beta: f64, y: Vector, z: Vector) void {
    cblasfeo.blasfeo_daxpby(x.len, alpha, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
}
pub fn vecmul(x: Vector, y: Vector, z: Vector) void {
    cblasfeo.blasfeo_dvecmul(x.len, x.dvec, x.start, y.dvec, y.start, z.dvec, z.start);
}
pub fn vecmulacc(x: Vector, y: Vector, z: Vector) void {
    cblasfeo.blasfeo_dvecmulacc(x.len, x.dvec, x.start, y.dvec, y.start, z.dvec, z.start);
}
pub fn vecmuldot(x: Vector, y: Vector, z: Vector) double {
    cblasfeo.blasfeo_dvecmuldot(x.len, x.dvec, x.start, y.dvec, y.start, z.dvec, z.start);
}
pub fn dot(x: Vector, y: Vector) double {
    cblasfeo.blasfeo_ddot(x.len, x.dvec, x.start, y.dvec, y.start);
}
pub fn rotg(a: f64, b: f64, c: *f64, s: *f64) void {
    cblasfeo.blasfeo_drotg(a, b, c, s);
}
pub fn colrot(A: Matrix, aj0: c_int, aj1: c_int, c: f64, s: f64) void {
    cblasfeo.blasfeo_dcolrot(A.rows, A.dmat, A.row_start, aj0, aj1, c, s);
}
pub fn rowrot(A: Matrix, ai0: c_int, ai1: c_int, c: f64, s: f64) void {
    cblasfeo.blasfeo_drowrot(A.rows, A.dmat, ai0, ai1, A.col_start, c, s);
}
pub fn gemv_n(alpha: f64, A: Matrix, x: Vector, beta: f64, y: Vector, z: Vector) void {
    cblasfeo.blasfeo_dgemv_n(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
}
pub fn gemv_t(alpha: f64, A: Matrix, x: Vector, beta: f64, y: Vector, z: Vector) void {
    cblasfeo.blasfeo_dgemv_t(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
}
pub fn trsv_lnn_mn(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrsv_lnn_mn(A.rows, A.cols, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trsv_ltn_mn(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrsv_ltn_mn(A.rows, A.cols, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trsv_lnn(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrsv_lnn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trsv_lnu(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrsv_lnu(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trsv_ltn(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrsv_ltn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trsv_ltu(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrsv_ltu(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trsv_unn(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrsv_unn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trsv_utn(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrsv_utn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trmv_lnn(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrmv_lnn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trmv_lnu(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrmv_lnu(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trmv_ltn(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrmv_ltn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trmv_ltu(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrmv_ltu(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trmv_unn(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrmv_unn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn trmv_utn(A: Matrix, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dtrmv_utn(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start, z.dvec, z.start);
}
pub fn gemv_nt(alpha_n: f64, alpha_t: f64, A: Matrix, x_n: Vector, x_t: Vector, beta_n: f64, beta_t: f64, y_n: Vector, y_t: Vector, z_n: Vector, z_t: Vector) void {
    cblasfeo.blasfeo_dgemv_nt(A.rows, A.cols, alpha_n, alpha_t, A.dmat, A.row_start, A.col_start, x_n.dvec, x_n.start, x_t.dvec, x_t.start, beta_n, beta_t, y_n.dvec, y_n.start, y_t.dvec, y_t.start, z_n.dvec, z_n.start, z_t.dvec, z_t.start);
}
pub fn symv_l(alpha: f64, A: Matrix, x: Vector, beta: f64, y: Vector, z: Vector) void {
    cblasfeo.blasfeo_dsymv_l(A.rows, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
}
pub fn symv_l_mn(alpha: f64, A: Matrix, x: Vector, beta: f64, y: Vector, z: Vector) void {
    cblasfeo.blasfeo_dsymv_l_mn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
}
pub fn symv_u(alpha: f64, A: Matrix, x: Vector, beta: f64, y: Vector, z: Vector) void {
    cblasfeo.blasfeo_dsymv_u(A.rows, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
}
pub fn ger(alpha: f64, x: Vector, y: Vector, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dger(C.rows, C.cols, alpha, x.dvec, x.start, y.dvec, y.start, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn gemv_d(alpha: f64, A: Vector, ai: c_int, x: Vector, beta: f64, y: Vector, z: Vector) void {
    cblasfeo.blasfeo_dgemv_d(A.len, alpha, A.dvec, ai, x.dvec, x.start, beta, y.dvec, y.start, z.dvec, z.start);
}
pub fn gemm_nn(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dgemm_nn(A.rows, A.cols, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn gemm_nt(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dgemm_nt(A.rows, A.cols, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn gemm_tn(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dgemm_tn(A.rows, A.cols, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn gemm_tt(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dgemm_tt(A.rows, A.cols, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn syrk_ln(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dsyrk_ln(A.rows, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn syrk_ln_mn(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dsyrk_ln_mn(A.rows, A.cols, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn syrk_lt(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dsyrk_lt(A.rows, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn syrk_un(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dsyrk_un(A.rows, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn syrk_ut(alpha: f64, A: Matrix, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dsyrk_ut(A.rows, D.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_llnn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_llnn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_llnu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_llnu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_lltn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_lltn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_lltu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_lltu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_lunn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_lunn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_lunu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_lunu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_lutn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_lutn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_lutu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_lutu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_rlnn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_rlnn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_rlnu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_rlnu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_rltn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_rltn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_rltu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_rltu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_runn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_runn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_runu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_runu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_rutn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_rutn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trmm_rutu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrmm_rutu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_llnn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_llnn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_llnu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_llnu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_lltn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_lltn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_lltu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_lltu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_lunn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_lunn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_lunu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_lunu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_lutn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_lutn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_lutu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_lutu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_rlnn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_rlnn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_rlnu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_rlnu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_rltn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_rltn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_rltu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_rltu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_runn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_runn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_runu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_runu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_rutn(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_rutn(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn trsm_rutu(alpha: f64, A: Matrix, B: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dtrsm_rutu(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn gemm_dn(alpha: f64, A: Vector, ai: c_int, B: Matrix, beta: f64, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dgemm_dn(B.rows, B.cols, alpha, A.dvec, ai, B.dmat, B.row_start, B.col_start, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn gemm_nd(alpha: f64, A: Matrix, B: Vector, bi: c_int, beta: f64, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dgemm_nd(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dvec, bi, beta, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn potrf_l(C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dpotrf_l(C.rows, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn potrf_l_mn(C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dpotrf_l_mn(C.rows, C.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn potrf_u(C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dpotrf_u(C.rows, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn syrk_dpotrf_ln(A: Matrix, B: Matrix, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dsyrk_dpotrf_ln(A.rows, D.cols, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn syrk_dpotrf_ln_mn(A: Matrix, B: Matrix, C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dsyrk_dpotrf_ln_mn(A.rows, A.cols, D.cols, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn getrf_np(C: Matrix, D: Matrix) void {
    cblasfeo.blasfeo_dgetrf_np(C.rows, C.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start);
}
pub fn getrf_rp(C: Matrix, D: Matrix, ipiv: *c_int) void {
    cblasfeo.blasfeo_dgetrf_rp(C.rows, C.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start, ipiv);
}
pub fn geqrf(C: Matrix, D: Matrix, work: *anyopaque) void {
    cblasfeo.blasfeo_dgeqrf(C.rows, C.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start, work);
}
pub fn orglq(C: Matrix, D: Matrix, work: *anyopaque) void {
    cblasfeo.blasfeo_dorglq(C.rows, C.cols, D.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start, work);
}
pub fn gelqf(C: Matrix, D: Matrix, work: *anyopaque) void {
    cblasfeo.blasfeo_dgelqf(C.rows, C.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start, work);
}
pub fn gelqf_pd(C: Matrix, D: Matrix, work: *anyopaque) void {
    cblasfeo.blasfeo_dgelqf_pd(C.rows, C.cols, C.dmat, C.row_start, C.col_start, D.dmat, D.row_start, D.col_start, work);
}
pub fn gelqf_pd_la(n1: c_int, L: Matrix, A: Matrix, work: *anyopaque) void {
    cblasfeo.blasfeo_dgelqf_pd_la(L.rows, n1, L.dmat, L.row_start, L.col_start, A.dmat, A.row_start, A.col_start, work);
}
pub fn gelqf_pd_lla(n1: c_int, L0: Matrix, L1: Matrix, A: Matrix, work: *anyopaque) void {
    cblasfeo.blasfeo_dgelqf_pd_lla(L0.rows, n1, L0.dmat, L0.row_start, L0.col_start, L1.dmat, L1.row_start, L1.col_start, A.dmat, A.row_start, A.col_start, work);
}

pub fn create_dmat(A: Matrix, memory: *anyopaque) void {
    cblasfeo.blasfeo_create_dmat(A.rows, A.cols, A.dmat, memory);
}
pub fn create_dvec(A: Vector, memory: *anyopaque) void {
    cblasfeo.blasfeo_create_dvec(A.len, A.dvec, memory);
}
pub fn pack_dmat(A: *f64, lda: c_int, B: Matrix) void {
    cblasfeo.blasfeo_pack_dmat(B.rows, B.cols, A, lda, B.dmat, B.row_start, B.col_start);
}
pub fn pack_l_dmat(A: *f64, lda: c_int, B: Matrix) void {
    cblasfeo.blasfeo_pack_l_dmat(B.rows, B.cols, A, lda, B.dmat, B.row_start, B.col_start);
}
pub fn pack_u_dmat(A: *f64, lda: c_int, B: Matrix) void {
    cblasfeo.blasfeo_pack_u_dmat(B.rows, B.cols, A, lda, B.dmat, B.row_start, B.col_start);
}
pub fn pack_tran_dmat(A: *f64, lda: c_int, B: Matrix) void {
    cblasfeo.blasfeo_pack_tran_dmat(B.rows, B.cols, A, lda, B.dmat, B.row_start, B.col_start);
}
pub fn pack_dvec(x: *f64, xi: c_int, y: Vector) void {
    cblasfeo.blasfeo_pack_dvec(y.len, x, xi, y.dvec, y.start);
}
pub fn unpack_dmat(A: Matrix, B: *f64, ldb: c_int) void {
    cblasfeo.blasfeo_unpack_dmat(A.rows, A.cols, A.dmat, A.row_start, A.col_start, B, ldb);
}
pub fn unpack_tran_dmat(A: Matrix, B: *f64, ldb: c_int) void {
    cblasfeo.blasfeo_unpack_tran_dmat(A.rows, A.cols, A.dmat, A.row_start, A.col_start, B, ldb);
}
pub fn unpack_dvec(x: Vector, y: *f64, yi: c_int) void {
    cblasfeo.blasfeo_unpack_dvec(x.len, x.dvec, x.start, y, yi);
}
pub fn dgese(alpha: f64, A: Matrix) void {
    cblasfeo.blasfeo_dgese(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start);
}
pub fn dgecp(A: Matrix, B: Matrix) void {
    cblasfeo.blasfeo_dgecp(A.rows, A.cols, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
}
pub fn dgesc(alpha: f64, A: Matrix) void {
    cblasfeo.blasfeo_dgesc(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start);
}
pub fn dgecpsc(alpha: f64, A: Matrix, B: Matrix) void {
    cblasfeo.blasfeo_dgecpsc(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
}
pub fn dtrcp_l(A: Matrix, B: Matrix) void {
    cblasfeo.blasfeo_dtrcp_l(A.rows, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
}
pub fn dtrcpsc_l(alpha: f64, A: Matrix, B: Matrix) void {
    cblasfeo.blasfeo_dtrcpsc_l(A.rows, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
}
pub fn dtrsc_l(alpha: f64, A: Matrix) void {
    cblasfeo.blasfeo_dtrsc_l(A.rows, alpha, A.dmat, A.row_start, A.col_start);
}
pub fn dgead(alpha: f64, A: Matrix, B: Matrix) void {
    cblasfeo.blasfeo_dgead(A.rows, A.cols, alpha, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
}
pub fn dvecad(alpha: f64, x: Vector, y: Vector) void {
    cblasfeo.blasfeo_dvecad(x.len, alpha, x.dvec, x.start, y.dvec, y.start);
}
pub fn dgetr(A: Matrix, B: Matrix) void {
    cblasfeo.blasfeo_dgetr(A.rows, A.cols, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
}
pub fn dtrtr_l(A: Matrix, B: Matrix) void {
    cblasfeo.blasfeo_dtrtr_l(A.rows, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
}
pub fn dtrtr_u(A: Matrix, B: Matrix) void {
    cblasfeo.blasfeo_dtrtr_u(A.rows, A.dmat, A.row_start, A.col_start, B.dmat, B.row_start, B.col_start);
}
pub fn ddiare(alpha: f64, A: Matrix) void {
    cblasfeo.blasfeo_ddiare(A.rows, alpha, A.dmat, A.row_start, A.col_start);
}
pub fn ddiain(alpha: f64, x: Vector, A: Matrix) void {
    cblasfeo.blasfeo_ddiain(A.rows, alpha, x.dvec, x.start, A.dmat, A.row_start, A.col_start);
}
pub fn ddiain_sp(alpha: f64, x: Vector, idx: *c_int, D: Matrix) void {
    cblasfeo.blasfeo_ddiain_sp(D.rows, alpha, x.dvec, x.start, idx, D.dmat, D.row_start, D.col_start);
}
pub fn ddiaex(alpha: f64, A: Matrix, x: Vector) void {
    cblasfeo.blasfeo_ddiaex(A.rows, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start);
}
pub fn ddiaex_sp(alpha: f64, idx: *c_int, D: Matrix, x: Vector) void {
    cblasfeo.blasfeo_ddiaex_sp(D.rows, alpha, idx, D.dmat, D.row_start, D.col_start, x.dvec, x.start);
}
pub fn ddiaad(alpha: f64, x: Vector, A: Matrix) void {
    cblasfeo.blasfeo_ddiaad(A.rows, alpha, x.dvec, x.start, A.dmat, A.row_start, A.col_start);
}
pub fn ddiaad_sp(alpha: f64, x: Vector, idx: *c_int, D: Matrix) void {
    cblasfeo.blasfeo_ddiaad_sp(D.rows, alpha, x.dvec, x.start, idx, D.dmat, D.row_start, D.col_start);
}
pub fn ddiaadin_sp(alpha: f64, x: Vector, y: Vector, idx: *c_int, D: Matrix) void {
    cblasfeo.blasfeo_ddiaadin_sp(D.rows, alpha, x.dvec, x.start, y.dvec, y.start, idx, D.dmat, D.row_start, D.col_start);
}
pub fn drowin(alpha: f64, x: Vector, A: Matrix) void {
    cblasfeo.blasfeo_drowin(A.rows, alpha, x.dvec, x.start, A.dmat, A.row_start, A.col_start);
}
pub fn drowex(alpha: f64, A: Matrix, x: Vector) void {
    cblasfeo.blasfeo_drowex(A.rows, alpha, A.dmat, A.row_start, A.col_start, x.dvec, x.start);
}
pub fn drowad(alpha: f64, x: Vector, A: Matrix) void {
    cblasfeo.blasfeo_drowad(A.rows, alpha, x.dvec, x.start, A.dmat, A.row_start, A.col_start);
}
pub fn drowad_sp(alpha: f64, x: Vector, idx: *c_int, D: Matrix) void {
    cblasfeo.blasfeo_drowad_sp(D.rows, alpha, x.dvec, x.start, idx, D.dmat, D.row_start, D.col_start);
}
pub fn drowsw(A: Matrix, C: Matrix) void {
    cblasfeo.blasfeo_drowsw(A.rows, A.dmat, A.row_start, A.col_start, C.dmat, C.row_start, C.col_start);
}
pub fn drowpe(ipiv: *c_int, A: Matrix) void {
    cblasfeo.blasfeo_drowpe(A.rows, ipiv, A.dmat);
}
pub fn drowpei(ipiv: *c_int, A: Matrix) void {
    cblasfeo.blasfeo_drowpei(A.rows, ipiv, A.dmat);
}
pub fn dcolex(A: Matrix, x: Vector) void {
    cblasfeo.blasfeo_dcolex(A.rows, A.dmat, A.row_start, A.col_start, x.dvec, x.start);
}
pub fn dcolin(x: Vector, A: Matrix) void {
    cblasfeo.blasfeo_dcolin(A.rows, x.dvec, x.start, A.dmat, A.row_start, A.col_start);
}
pub fn dcolad(alpha: f64, x: Vector, A: Matrix) void {
    cblasfeo.blasfeo_dcolad(A.rows, alpha, x.dvec, x.start, A.dmat, A.row_start, A.col_start);
}
pub fn dcolsc(alpha: f64, A: Matrix) void {
    cblasfeo.blasfeo_dcolsc(A.rows, alpha, A.dmat, A.row_start, A.col_start);
}
pub fn dcolsw(A: Matrix, C: Matrix) void {
    cblasfeo.blasfeo_dcolsw(A.rows, A.dmat, A.row_start, A.col_start, C.dmat, C.row_start, C.col_start);
}
pub fn dcolpe(ipiv: *c_int, A: Matrix) void {
    cblasfeo.blasfeo_dcolpe(A.rows, ipiv, A.dmat);
}
pub fn dcolpei(ipiv: *c_int, A: Matrix) void {
    cblasfeo.blasfeo_dcolpei(A.rows, ipiv, A.dmat);
}
pub fn dvecse(alpha: f64, x: Vector) void {
    cblasfeo.blasfeo_dvecse(x.len, alpha, x.dvec, x.start);
}
pub fn dveccp(x: Vector, y: Vector) void {
    cblasfeo.blasfeo_dveccp(x.len, x.dvec, x.start, y.dvec, y.start);
}
pub fn dvecsc(alpha: f64, x: Vector) void {
    cblasfeo.blasfeo_dvecsc(x.len, alpha, x.dvec, x.start);
}
pub fn dveccpsc(alpha: f64, x: Vector, y: Vector) void {
    cblasfeo.blasfeo_dveccpsc(x.len, alpha, x.dvec, x.start, y.dvec, y.start);
}
pub fn dvecad_sp(alpha: f64, x: Vector, idx: *c_int, z: Vector) void {
    cblasfeo.blasfeo_dvecad_sp(x.len, alpha, x.dvec, x.start, idx, z.dvec, z.start);
}
pub fn dvecin_sp(alpha: f64, x: Vector, idx: *c_int, z: Vector) void {
    cblasfeo.blasfeo_dvecin_sp(x.len, alpha, x.dvec, x.start, idx, z.dvec, z.start);
}
pub fn dvecex_sp(alpha: f64, idx: *c_int, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dvecex_sp(x.len, alpha, idx, x.dvec, x.start, z.dvec, z.start);
}
pub fn dvecexad_sp(alpha: f64, idx: *c_int, x: Vector, z: Vector) void {
    cblasfeo.blasfeo_dvecexad_sp(x.len, alpha, idx, x.dvec, x.start, z.dvec, z.start);
}
pub fn dvecze(m: Vector, v: Vector, e: Vector) void {
    cblasfeo.blasfeo_dvecze(m.len, m.dvec, m.start, v.dvec, v.start, e.dvec, e.start);
}
pub fn dvecnrm_inf(x: Vector, ptr_norm: *f64) void {
    cblasfeo.blasfeo_dvecnrm_inf(x.len, x.dvec, x.start, ptr_norm);
}
pub fn dvecpe(ipiv: *c_int, x: Vector) void {
    cblasfeo.blasfeo_dvecpe(x.len, ipiv, x.dvec, x.start);
}
pub fn dvecpei(ipiv: *c_int, x: Vector) void {
    cblasfeo.blasfeo_dvecpei(x.len, ipiv, x.dvec, x.start);
}
