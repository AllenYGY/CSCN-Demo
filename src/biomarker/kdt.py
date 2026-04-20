import copy
import random


def qnth_element(a, l, r, k, dim):
    if l == r:
        return a[l]
    pivot_idx = random.randint(l, r)
    x = a[pivot_idx][dim]
    i = l - 1
    j = r + 1
    while i < j:
        i += 1
        j -= 1
        while i <= r and a[i][dim] < x:
            i += 1
        while j >= l and a[j][dim] > x:
            j -= 1
        if i < j:
            a[i], a[j] = a[j], a[i]
    if k <= j:
        return qnth_element(a, l, j, k, dim)
    return qnth_element(a, j + 1, r, k, dim)


class KDT_Node:
    __slots__ = ("point", "left", "right", "low_bounds", "up_bounds", "size")

    def __init__(self, point, left, right):
        self.point = point
        self.left = left
        self.right = right

    def push_up(self):
        self.low_bounds = self.point.copy()
        self.up_bounds = self.point.copy()
        self.size = 1
        for child in (self.left, self.right):
            if child is not None:
                self.size += child.size
                for i in range(len(self.point)):
                    self.low_bounds[i] = min(self.low_bounds[i], child.low_bounds[i])
                    self.up_bounds[i] = max(self.up_bounds[i], child.up_bounds[i])


class KDT:
    __slots__ = ("points", "K", "root")

    def __init__(self, points):
        self.points = copy.deepcopy(points)
        self.K = len(points[0])
        self.root = self.build(0, len(points) - 1, 0)

    def build(self, l, r, axis):
        if l > r:
            return None
        axis %= self.K
        n = r - l + 1
        idx = n // 2 + l
        point = qnth_element(self.points, l, r, idx, axis)
        rgt = self.build(idx + 1, r, axis + 1)
        lft = self.build(l, idx - 1, axis + 1)
        node = KDT_Node(point, lft, rgt)
        node.push_up()
        return node

    def query(self, node, q):
        if not node:
            return 0
        flag = False
        for i in range(self.K):
            low = q[i][0]
            up = q[i][1]
            if up < low:
                continue
            if not (low <= node.low_bounds[i] and up >= node.up_bounds[i]):
                flag = True
                break
        if not flag:
            return node.size
        for i in range(self.K):
            low = q[i][0]
            up = q[i][1]
            if up < low:
                continue
            if node.low_bounds[i] > up or node.up_bounds[i] < low:
                return 0
        ans = 0
        flag = False
        for i in range(self.K):
            low = q[i][0]
            up = q[i][1]
            if up < low:
                continue
            if not (node.point[i] <= up and node.point[i] >= low):
                flag = True
                break
        if not flag:
            ans += 1
        ans += self.query(node.left, q)
        ans += self.query(node.right, q)
        return ans

    def query_cnt(self, q):
        return self.query(self.root, q)
