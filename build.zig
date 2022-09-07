const std = @import("std");

pub fn dependOnBlasfeo(step: *std.build.LibExeObjStep) void {
    // add blasfeo system dependency
    step.addIncludeDir("/opt/blasfeo/include");
    step.addLibPath("/opt/blasfeo/lib");
    step.linkSystemLibrary(":libblasfeo.a");
    step.linkLibC();
}

pub fn build(b: *std.build.Builder) void {
    // Standard target options allows the person running `zig build` to choose
    // what target to build for. Here we do not override the defaults, which
    // means any target is allowed, and the default is native. Other options
    // for restricting supported target set are available.
    const target = b.standardTargetOptions(.{});

    // Standard release options allow the person running `zig build` to select
    // between Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall.
    const mode = b.standardReleaseOptions();

    const blasfeo_tests = b.addTest("src/blasfeo.zig");
    blasfeo_tests.setTarget(target);
    blasfeo_tests.setBuildMode(mode);
    dependOnBlasfeo(blasfeo_tests);
    const blasfeo_tests_step = b.step("test-blasfeo", "Run blasfeo tests");
    blasfeo_tests_step.dependOn(&blasfeo_tests.step);

    const factor_graph_tests = b.addTest("src/factor_graph.zig");
    factor_graph_tests.setTarget(target);
    factor_graph_tests.setBuildMode(mode);
    dependOnBlasfeo(factor_graph_tests);
    const factor_graph_tests_step = b.step("test-factor-graph", "Run blasfeo tests");
    factor_graph_tests_step.dependOn(&factor_graph_tests.step);

    const upriority_queue_tests = b.addTest("src/upriority_queue.zig");
    upriority_queue_tests.setTarget(target);
    upriority_queue_tests.setBuildMode(mode);
    dependOnBlasfeo(upriority_queue_tests);
    const upriority_queue_tests_step = b.step("test-upq", "Run updatable priority queue");
    upriority_queue_tests_step.dependOn(&upriority_queue_tests.step);
}
