const std = @import("std");
const Allocator = std.mem.Allocator;
const assert = std.debug.assert;
const Order = std.math.Order;
const testing = std.testing;
const expect = testing.expect;
const expectEqual = testing.expectEqual;
const expectError = testing.expectError;

pub const Priority = struct {
    id: usize,
    priority: f64,
};

fn compareFn(a: Priority, b: Priority) Order {
    if (a.priority < b.priority) return Order.lt;
    if (a.priority > b.priority) return Order.gt;
    return Order.eq;
}

/// Priority queue for storing ids associated with a constant time
/// updateable priority. Based heavily on std.PriorityQueue
pub const UPriorityQueue = struct {
    const Self = @This();

    items: []Priority,
    len: usize,
    allocator: Allocator,
    idx_lookup: std.AutoHashMap(usize, usize),

    /// Initialize and return a priority queue.
    pub fn init(allocator: Allocator) Self {
        // zig fmt: off
        return Self{
            .items = &[_]Priority{},
            .len = 0,
            .allocator = allocator,
            .idx_lookup = std.AutoHashMap(usize, usize).init(allocator)
        };
        // zig fmt: on
    }

    /// Free memory used by the queue.
    pub fn deinit(self: *Self) void {
        self.allocator.free(self.items);
        self.idx_lookup.deinit();
    }

    /// Insert a new element, maintaining priority.
    pub fn add(self: *Self, id: usize, priority: f64) !void {
        try self.ensureUnusedCapacity(1);
        addUnchecked(self, id, priority);
    }

    fn addUnchecked(self: *Self, id: usize, priority: f64) void {
        self.items[self.len] = Priority{ .id = id, .priority = priority };
        self.idx_lookup.putAssumeCapacity(id, self.len);
        siftUp(self, self.len);
        self.len += 1;
    }

    fn siftUp(self: *Self, start_index: usize) void {
        var child_index = start_index;
        while (child_index > 0) {
            var parent_index = ((child_index - 1) >> 1);
            const child = self.items[child_index];
            const parent = self.items[parent_index];

            if (compareFn(child, parent) != .lt) break;

            self.items[parent_index] = child;
            self.items[child_index] = parent;

            self.idx_lookup.put(child.id, parent_index) catch unreachable;
            self.idx_lookup.put(parent.id, child_index) catch unreachable;

            child_index = parent_index;
        }
    }

    // /// Add each element in `items` to the queue.
    // pub fn addSlice(self: *Self, items: []const Priority) !void {
    //     try self.ensureUnusedCapacity(items.len);
    //     for (items) |e| {
    //         self.addUnchecked(e);
    //     }
    // }

    /// Look at the highest priority element in the queue. Returns
    /// `null` if empty.
    pub fn peek(self: *Self) ?Priority {
        return if (self.len > 0) self.items[0] else null;
    }

    /// Pop the highest priority element from the queue. Returns
    /// `null` if empty.
    pub fn removeOrNull(self: *Self) ?Priority {
        return if (self.len > 0) self.remove() else null;
    }

    /// Remove and return the highest priority element from the
    /// queue.
    pub fn remove(self: *Self) Priority {
        return self.removeIndex(0);
    }

    /// Remove and return element at index. Indices are in the
    /// same order as iterator, which is not necessarily priority
    /// order.
    pub fn removeIndex(self: *Self, index: usize) Priority {
        assert(self.len > index);
        const last = self.items[self.len - 1];
        const item = self.items[index];
        self.items[index] = last;
        self.len -= 1;

        if (index == 0) {
            siftDown(self, index);
        } else {
            const parent_index = ((index - 1) >> 1);
            const parent = self.items[parent_index];
            if (compareFn(last, parent) == .gt) {
                siftDown(self, index);
            } else {
                siftUp(self, index);
            }
        }

        return item;
    }

    /// Return the number of elements remaining in the priority
    /// queue.
    pub fn count(self: Self) usize {
        return self.len;
    }

    /// Return the number of elements that can be added to the
    /// queue before more memory is allocated.
    pub fn capacity(self: Self) usize {
        return self.items.len;
    }

    fn siftDown(self: *Self, start_index: usize) void {
        var index = start_index;
        const half = self.len >> 1;
        while (true) {
            var left_index = (index << 1) + 1;
            var right_index = left_index + 1;
            var left = if (left_index < self.len) self.items[left_index] else null;
            var right = if (right_index < self.len) self.items[right_index] else null;

            var smallest_index = index;
            var smallest = self.items[index];

            if (left) |e| {
                if (compareFn(e, smallest) == .lt) {
                    smallest_index = left_index;
                    smallest = e;
                }
            }

            if (right) |e| {
                if (compareFn(e, smallest) == .lt) {
                    smallest_index = right_index;
                    smallest = e;
                }
            }

            if (smallest_index == index) return;

            self.items[smallest_index] = self.items[index];
            self.items[index] = smallest;

            self.idx_lookup.put(self.items[smallest_index].id, smallest_index) catch unreachable;
            self.idx_lookup.put(self.items[index].id, index) catch unreachable;

            index = smallest_index;

            if (index >= half) return;
        }
    }

    // /// PriorityQueue takes ownership of the passed in slice. The slice must have been
    // /// allocated with `allocator`.
    // /// Deinitialize with `deinit`.
    // pub fn fromOwnedSlice(allocator: Allocator, items: []Priority) Self {
    //     // zig fmt: off
    //     var queue = Self{
    //         .items = items,
    //         .len = items.len,
    //         .allocator = allocator,
    //         .idx_lookup = std.AutoHashMap(Priority, usize).init(allocator)
    //     };
    //     // zig fmt: on

    //     if (queue.len <= 1) return queue;

    //     const half = (queue.len >> 1) - 1;
    //     var i: usize = 0;
    //     while (i <= half) : (i += 1) {
    //         queue.siftDown(half - i);
    //     }
    //     return queue;
    // }

    pub const ensureCapacity = @compileError("deprecated; use ensureUnusedCapacity or ensureTotalCapacity");

    /// Ensure that the queue can fit at least `new_capacity` items.
    pub fn ensureTotalCapacity(self: *Self, new_capacity: usize) !void {
        var better_capacity = self.capacity();
        if (better_capacity >= new_capacity) return;
        while (true) {
            better_capacity += better_capacity / 2 + 8;
            if (better_capacity >= new_capacity) break;
        }
        self.items = try self.allocator.realloc(self.items, better_capacity);
        try self.idx_lookup.ensureTotalCapacity(@intCast(u32, better_capacity));
    }

    /// Ensure that the queue can fit at least `additional_count` **more** item.
    pub fn ensureUnusedCapacity(self: *Self, additional_count: usize) !void {
        return self.ensureTotalCapacity(self.len + additional_count);
    }

    /// Reduce allocated capacity to `new_len`.
    pub fn shrinkAndFree(self: *Self, new_len: usize) void {
        assert(new_len <= self.items.len);

        // Cannot shrink to smaller than the current queue size without invalidating the heap property
        assert(new_len >= self.len);

        self.items = self.allocator.realloc(self.items[0..], new_len) catch |e| switch (e) {
            error.OutOfMemory => { // no problem, capacity is still correct then.
                self.items.len = new_len;
                return;
            },
        };
    }

    pub fn update(self: *Self, id: usize, new_priority: f64) !void {
        const update_index = self.idx_lookup.get(id) orelse return error.ElementNotFound;
        const old_elem: Priority = self.items[update_index];
        self.items[update_index].priority = new_priority;
        const new_elem: Priority = self.items[update_index];
        switch (compareFn(new_elem, old_elem)) {
            .lt => siftUp(self, update_index),
            .gt => siftDown(self, update_index),
            .eq => {}, // Nothing to do as the items have equal priority
        }
    }

    pub const Iterator = struct {
        queue: *UPriorityQueue,
        count: usize,

        pub fn next(it: *Iterator) ?Priority {
            if (it.count >= it.queue.len) return null;
            const out = it.count;
            it.count += 1;
            return it.queue.items[out];
        }

        pub fn reset(it: *Iterator) void {
            it.count = 0;
        }
    };

    /// Return an iterator that walks the queue without consuming
    /// it. Invalidated if the heap is modified.
    pub fn iterator(self: *Self) Iterator {
        return Iterator{
            .queue = self,
            .count = 0,
        };
    }

    fn dump(self: *Self) void {
        const print = std.debug.print;
        print("{{ ", .{});
        print("items: ", .{});
        for (self.items) |e, i| {
            if (i >= self.len) break;
            print("{}, ", .{e});
        }
        print("array: ", .{});
        for (self.items) |e| {
            print("{}, ", .{e});
        }
        print("len: {} ", .{self.len});
        print("capacity: {}", .{self.capacity()});
        print(" }}\n", .{});
    }
};

const PQlt = UPriorityQueue;

test "upriorityQueue: add and remove min heap" {
    var queue = PQlt.init(testing.allocator);
    defer queue.deinit();

    try queue.add(54, 54);
    try queue.add(12, 12);
    try queue.add(7, 7);
    try queue.add(23, 23);
    try queue.add(25, 25);
    try queue.add(13, 13);
    try std.testing.expect(queue.remove().id == 7);
    try std.testing.expect(queue.remove().id == 12);
    try std.testing.expect(queue.remove().id == 13);
    try std.testing.expect(queue.remove().id == 23);
    try std.testing.expect(queue.remove().id == 25);
    try std.testing.expect(queue.remove().id == 54);
}

test "upriorityQueue: add and remove same min heap" {
    var queue = PQlt.init(testing.allocator);
    defer queue.deinit();

    try queue.add(1, 1);
    try queue.add(1, 1);
    try queue.add(2, 2);
    try queue.add(2, 2);
    try queue.add(1, 1);
    try queue.add(1, 1);
    try std.testing.expect(queue.remove().id == 1);
    try std.testing.expect(queue.remove().id == 1);
    try std.testing.expect(queue.remove().id == 1);
    try std.testing.expect(queue.remove().id == 1);
    try std.testing.expect(queue.remove().id == 2);
    try std.testing.expect(queue.remove().id == 2);
}

test "upriorityQueue: removeOrNull on empty" {
    var queue = PQlt.init(testing.allocator);
    defer queue.deinit();

    try expect(queue.removeOrNull() == null);
}

test "upriorityQueue: edge case 3 elements" {
    var queue = PQlt.init(testing.allocator);
    defer queue.deinit();

    try queue.add(9, 9);
    try queue.add(3, 3);
    try queue.add(2, 2);
    try std.testing.expect(queue.remove().id == 2);
    try std.testing.expect(queue.remove().id == 3);
    try std.testing.expect(queue.remove().id == 9);
}

test "upriorityQueue: peek" {
    var queue = PQlt.init(testing.allocator);
    defer queue.deinit();

    try expect(queue.peek() == null);
    try queue.add(9, 9);
    try queue.add(3, 3);
    try queue.add(2, 2);
    try std.testing.expect(queue.peek().?.id == 2);
    try std.testing.expect(queue.peek().?.id == 2);
}

test "upriorityQueue: sift up with odd indices" {
    var queue = PQlt.init(testing.allocator);
    defer queue.deinit();
    const items = [_]u32{ 15, 7, 21, 14, 13, 22, 12, 6, 7, 25, 5, 24, 11, 16, 15, 24, 2, 1 };
    for (items) |e| {
        try queue.add(e, @intToFloat(f64, e));
    }

    const sorted_items = [_]u32{ 1, 2, 5, 6, 7, 7, 11, 12, 13, 14, 15, 15, 16, 21, 22, 24, 24, 25 };
    for (sorted_items) |e| {
        try std.testing.expect(queue.remove().id == e);
    }
}

test "upriorityQueue: iterator" {
    var queue = PQlt.init(testing.allocator);
    var map = std.AutoHashMap(usize, void).init(testing.allocator);
    defer {
        queue.deinit();
        map.deinit();
    }

    const items = [_]u32{ 54, 12, 7, 23, 25, 13 };
    for (items) |e| {
        _ = try queue.add(e, @intToFloat(f64, e));
        try map.put(e, {});
    }

    var it = queue.iterator();
    while (it.next()) |e| {
        _ = map.remove(e.id);
    }

    try expectEqual(@as(usize, 0), map.count());
}

test "upriorityQueue: remove at index" {
    var queue = PQlt.init(testing.allocator);
    defer queue.deinit();

    const items = [_]u32{ 2, 1, 8, 9, 3, 4, 5 };
    for (items) |e| {
        _ = try queue.add(e, @intToFloat(f64, e));
    }

    var it = queue.iterator();
    var idx: usize = 0;
    const two_idx = while (it.next()) |elem| {
        if (elem.id == 2)
            break idx;
        idx += 1;
    } else unreachable;
    var sorted_items = [_]u32{ 1, 3, 4, 5, 8, 9 };
    try expectEqual(queue.removeIndex(two_idx).id, 2);

    var i: usize = 0;
    while (queue.removeOrNull()) |n| : (i += 1) {
        try expectEqual(n.id, sorted_items[i]);
    }
    try expectEqual(queue.removeOrNull(), null);
}

test "upriorityQueue: iterator while empty" {
    var queue = PQlt.init(testing.allocator);
    defer queue.deinit();

    var it = queue.iterator();

    try expectEqual(it.next(), null);
}

test "upriorityQueue: shrinkAndFree" {
    var queue = PQlt.init(testing.allocator);
    defer queue.deinit();

    try queue.ensureTotalCapacity(4);
    try expect(queue.capacity() >= 4);

    try queue.add(1, 1);
    try queue.add(2, 2);
    try queue.add(3, 3);
    try expect(queue.capacity() >= 4);
    try expectEqual(@as(usize, 3), queue.len);

    queue.shrinkAndFree(3);
    try expectEqual(@as(usize, 3), queue.capacity());
    try expectEqual(@as(usize, 3), queue.len);

    try std.testing.expect(queue.remove().id == 1);
    try std.testing.expect(queue.remove().id == 2);
    try std.testing.expect(queue.remove().id == 3);
    try expect(queue.removeOrNull() == null);
}

test "upriorityQueue: update min heap" {
    var queue = PQlt.init(testing.allocator);
    defer queue.deinit();

    try queue.add(55, 1);
    try queue.add(44, 2);
    try queue.add(11, 3);
    try queue.update(55, 3);
    try queue.update(44, 2);
    try queue.update(11, 1);

    try std.testing.expect(queue.remove().id == 11);
    try std.testing.expect(queue.remove().id == 44);
    try std.testing.expect(queue.remove().id == 55);
}

test "upriorityQueue: update same min heap" {
    var queue = PQlt.init(testing.allocator);
    defer queue.deinit();

    try queue.add(1, 10);
    try queue.add(2, 9);
    try queue.add(3, 3);
    try queue.add(4, 4);
    try queue.update(1, 1);
    try queue.update(2, 2);

    try std.testing.expect(queue.remove().id == 1);
    try std.testing.expect(queue.remove().id == 2);
    try std.testing.expect(queue.remove().id == 3);
    try std.testing.expect(queue.remove().id == 4);
}
