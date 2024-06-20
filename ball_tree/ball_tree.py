import numpy as np


class Node:
    def __init__(self, pivot, points, radius):
        self.pivot = pivot  # Центроид
        self.radius = radius
        self.points = points

    def __str__(self):
        self.stringList = ''.join(str(symbol) for symbol in self.points)
        return ("pivot = " + str(self.pivot) + " radius = " + str(self.radius) + " points = "
                + str(self.stringList))


class Point:
    def __init__(self, point, index):
        self.point = point
        self.index = index

    def __str__(self):
        return "Vector " + str(self.point) + " with index " + str(self.index)


class BallTree:
    def __init__(self, leaf_size=30):
        self.leaf_size = leaf_size
        self.nodes = {}
        self.kNSet = set()
        self.points = []

    def fit(self, x):
        self.points = []
        for i in range(len(x)):
            self.points.append(Point(x.iloc[i], i))

        self._construct_tree(self.leaf_size, 1, 0, len(self.points))
        return self

    def search(self, x, n_neighbors):
        D_for_all = []
        I_for_all = []
        for i in range(len(x)):
            self.kNSet.clear()
            self._search_ball_subtree(1, x.iloc[i], n_neighbors)
            D = []
            I = []

            for neighbour in self.kNSet:
                # append distance
                D.append(neighbour[0])
                # append index
                I.append(neighbour[1])
            D_for_all.append(D)
            I_for_all.append(I)
        return D_for_all, I_for_all

    def _construct_tree(self, leaf_size, vertex_index, start_point, end_point):
        if end_point - start_point <= 0:
            return

        dimension = self._get_max_spread_dimension(self.points[start_point:end_point])
        self.points = self.sorting_by_dimension(self.points, start_point, end_point, dimension)
        centroid_index = self._get_centroid_index(self.points[start_point:end_point], dimension) + start_point

        radius = self._get_max_distance(self.points[centroid_index], self.points[start_point:end_point])
        self.nodes[vertex_index] = Node(self.points[centroid_index], [], radius)
        if end_point - start_point <= leaf_size:
            for i in range(start_point, end_point):
                self.nodes[vertex_index].points.append(self.points[i])
            return
        else:
            self.nodes[vertex_index].points.append(self.points[centroid_index])
        self._construct_tree(leaf_size, vertex_index * 2, start_point, centroid_index)
        self._construct_tree(leaf_size, vertex_index * 2 + 1, centroid_index + 1, end_point)

    def _search_ball_subtree(self, vertex_index, new_point, n_neighbors):
        if len(self.nodes[vertex_index].points) > 1:
            for point in self.nodes[vertex_index].points:
                distance = self.distance(new_point, point.point)
                if len(self.kNSet) >= n_neighbors and max(self.kNSet)[0] > distance:
                    self.kNSet.remove(max(self.kNSet))
                if len(self.kNSet) < n_neighbors:
                    self.kNSet.add((distance, point.index))
            return
        left_child = vertex_index * 2
        right_child = vertex_index * 2 + 1
        distance = self.distance(new_point, self.nodes[vertex_index].pivot.point)

        if len(self.kNSet) >= n_neighbors and max(self.kNSet)[0] > distance:
            self.kNSet.remove(max(self.kNSet))
        if len(self.kNSet) < n_neighbors:
            self.kNSet.add((distance, self.nodes[vertex_index].pivot.index))

        self._check_tree_branch(left_child, new_point, n_neighbors)
        self._check_tree_branch(right_child, new_point, n_neighbors)

    def _get_max_distance(self, point, points):
        points = self.get_np_array_of_points(points)
        point = point.point
        # print(points)

        max = 0
        for i in range(len(points)):
            current_dist = self.distance(point, points[i])
            if current_dist > max:
                max = current_dist
        return max

    def _get_max_spread_dimension(self, x):
        x = self.get_np_array_of_points(x)
        difference = []
        for i in range(x.shape[1]):
            maximum = x[0][i]
            minimum = x[0][i]
            for j in range(x.shape[0]):
                maximum = max(maximum, x[j][i])
                minimum = min(minimum, x[j][i])
            difference.append((maximum - minimum, i))
        return max(difference)[1]

    def _get_centroid_index(self, x, dimension):
        x = self.get_np_array_of_points(x)
        if x.shape[0] == 0:
            print("Ошибка! Передан нулевой массив в поиск центроида")
        sum = 0
        for i in range(x.shape[0]):  # Размер x не нулевой
            sum = sum + x[i][dimension]
        sum = sum / x.shape[0]
        minimum = abs(sum - x[0][dimension])
        index = 0
        for i in range(x.shape[0]):
            if abs(x[i][dimension] - sum) < minimum:
                minimum = abs(x[i][dimension] - sum)
                index = i
        return index

    def _check_tree_branch(self, branch_index, new_point, n_neighbors):
        if branch_index in self.nodes and len(self.nodes[branch_index].points) != 0:
            distance = self.distance(new_point, self.nodes[branch_index].pivot.point)
            if len(self.kNSet) < n_neighbors or max(self.kNSet)[0] > distance - self.nodes[branch_index].radius:
                self._search_ball_subtree(branch_index, new_point, n_neighbors)

    @staticmethod
    def get_np_array_of_points(arr):
        x = []
        for i in range(len(arr)):
            x.append(arr[i].point)
        return np.array(x)

    @staticmethod
    def distance(a, b):
        if len(a) != len(b):
            print("Размерности объектов не одинаковые!")
            return "Error"
        return np.linalg.norm(a - b)

    @staticmethod
    def sorting_by_dimension(x, start_point, end_point, dimension):
        left = x[:start_point]
        right = x[end_point:]
        x = np.array(sorted(x[start_point:end_point], key=lambda point: point.point[dimension]))
        x = np.concatenate((left, x, right), 0)
        return x


# if __name__ == '__main__':
#     d = {}
#     d[56] = 1
#     d.
#     print(len(d[57]))