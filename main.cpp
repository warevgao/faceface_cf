#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <limits>
#include <set>
#include <tuple>
#include <algorithm>

// ==========================================
// 1. 基础几何库
// ==========================================
const double EPSILON = 1e-9;

class SPAposition {
private:
	double _x, _y, _z;

public:
	SPAposition(double x = 0, double y = 0, double z = 0) : _x(x), _y(y), _z(z) {}

	double x() const { return _x; }
	double y() const { return _y; }
	double z() const { return _z; }

	SPAposition operator+(const SPAposition& other) const { return SPAposition(_x + other._x, _y + other._y, _z + other._z); }
	SPAposition operator-(const SPAposition& other) const { return SPAposition(_x - other._x, _y - other._y, _z - other._z); }
	SPAposition operator*(double s) const { return SPAposition(_x * s, _y * s, _z * s); }
	SPAposition operator/(double s) const { return SPAposition(_x / s, _y / s, _z / s); }

	static double dot(const SPAposition& a, const SPAposition& b) {
		return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
	}

	static SPAposition cross(const SPAposition& a, const SPAposition& b) {
		return SPAposition(a.y() * b.z() - a.z() * b.y(), a.z() * b.x() - a.x() * b.z(), a.x() * b.y() - a.y() * b.x());
	}

	double magSq() const { return _x * _x + _y * _y + _z * _z; }
	double mag() const { return std::sqrt(magSq()); }

	SPAposition normalize() const {
		double m = mag();
		if (m < EPSILON) return *this;
		return *this / m;
	}
};

struct Triangle {
	SPAposition p0, p1, p2;
};

// ==========================================
// 2. 几何算法: 包含相交检测和最近点投影
// ==========================================
namespace Geometry {

	// 辅助：点到线段最近点
	SPAposition closestPointOnSegment(const SPAposition& p, const SPAposition& a, const SPAposition& b) {
		SPAposition ab = b - a;
		double t = SPAposition::dot(p - a, ab) / ab.magSq();
		t = std::max(0.0, std::min(1.0, t));
		return a + ab * t;
	}

	// 辅助：点到三角形最近点
	SPAposition closestPointOnTriangle(const SPAposition& p, const Triangle& tri) {
		SPAposition ab = tri.p1 - tri.p0;
		SPAposition ac = tri.p2 - tri.p0;
		SPAposition ap = p - tri.p0;

		// 简单法向量检查
		SPAposition n = SPAposition::cross(ab, ac);
		double nSq = n.magSq();
		SPAposition p_proj = p - n * (SPAposition::dot(ap, n) / nSq);

		// 检查重心坐标 (简化版: 检查是否在三边内侧)
		// 这里的简化版使用边缘平面测试
		SPAposition bc = tri.p2 - tri.p1;
		SPAposition ca = tri.p0 - tri.p2;

		if (SPAposition::dot(SPAposition::cross(ab, n), p_proj - tri.p0) > 0) goto check_edges;
		if (SPAposition::dot(SPAposition::cross(bc, n), p_proj - tri.p1) > 0) goto check_edges;
		if (SPAposition::dot(SPAposition::cross(ca, n), p_proj - tri.p2) > 0) goto check_edges;

		return p_proj; // 在内部

	check_edges:
		SPAposition c1 = closestPointOnSegment(p, tri.p0, tri.p1);
		SPAposition c2 = closestPointOnSegment(p, tri.p1, tri.p2);
		SPAposition c3 = closestPointOnSegment(p, tri.p2, tri.p0);

		double d1 = (p - c1).magSq();
		double d2 = (p - c2).magSq();
		double d3 = (p - c3).magSq();

		if (d1 <= d2 && d1 <= d3) return c1;
		if (d2 <= d1 && d2 <= d3) return c2;
		return c3;
	}

	// === 核心：三角面片求交 ===
	// 返回 true 表示相交，并填充 out1, out2 为交线段端点
	bool intersectTriTri(const Triangle& t1, const Triangle& t2, SPAposition& out1, SPAposition& out2) {
		// 1. 计算 T1 所在平面的方程 N1 * X + d1 = 0
		SPAposition n1 = SPAposition::cross(t1.p1 - t1.p0, t1.p2 - t1.p0).normalize();
		double d1 = -SPAposition::dot(n1, t1.p0);

		// 2. 计算 T2 顶点到 T1 平面的有符号距离
		double du0 = SPAposition::dot(n1, t2.p0) + d1;
		double du1 = SPAposition::dot(n1, t2.p1) + d1;
		double du2 = SPAposition::dot(n1, t2.p2) + d1;

		// 如果距离全为正或全为负，则不相交
		if (std::abs(du0) > EPSILON && std::abs(du1) > EPSILON && std::abs(du2) > EPSILON) {
			if ((du0 > 0 && du1 > 0 && du2 > 0) || (du0 < 0 && du1 < 0 && du2 < 0)) return false;
		}

		// 3. 计算 T2 所在平面方程 N2 * X + d2 = 0
		SPAposition n2 = SPAposition::cross(t2.p1 - t2.p0, t2.p2 - t2.p0).normalize();
		double d2 = -SPAposition::dot(n2, t2.p0);

		// 4. 计算 T1 顶点到 T2 平面的有符号距离
		double dv0 = SPAposition::dot(n2, t1.p0) + d2;
		double dv1 = SPAposition::dot(n2, t1.p1) + d2;
		double dv2 = SPAposition::dot(n2, t1.p2) + d2;

		if (std::abs(dv0) > EPSILON && std::abs(dv1) > EPSILON && std::abs(dv2) > EPSILON) {
			if ((dv0 > 0 && dv1 > 0 && dv2 > 0) || (dv0 < 0 && dv1 < 0 && dv2 < 0)) return false;
		}

		// 5. 计算交线的方向 D = N1 x N2
		SPAposition D = SPAposition::cross(n1, n2);

		// 如果 D 接近0，说明平面平行（共面情况暂忽略，视作不相交或由最近点处理）
		if (D.magSq() < EPSILON) return false;

		// 6. 核心逻辑：投影间隔 (Interval Overlap)
		// 将两个三角形沿 D 方向投影到一维直线上，求重叠区间
		// 这里简化实现：直接找 T1 穿过 Plane2 的线段，和 T2 穿过 Plane1 的线段

		auto get_inter_interval = [&](const Triangle& tri, double dist0, double dist1, double dist2, double& t_min, double& t_max) {
			// 找到跨越平面的两条边，计算 t 值
			// 投影公式: P = A + (B-A) * (-distA / (distB - distA))
			// 我们只需要在 Line 方向上的投影值

			// 这里的完整实现比较长，为了代码简洁，我们采用近似策略：
			// 找出 Tri 穿过平面的两点 P_start, P_end
			std::vector<SPAposition> pts;
			if (dist0 * dist1 <= 0) pts.push_back(tri.p0 + (tri.p1 - tri.p0) * (-dist0 / (dist1 - dist0)));
			if (dist1 * dist2 <= 0) pts.push_back(tri.p1 + (tri.p2 - tri.p1) * (-dist1 / (dist2 - dist1)));
			if (dist2 * dist0 <= 0) pts.push_back(tri.p2 + (tri.p0 - tri.p2) * (-dist2 / (dist0 - dist2)));

			if (pts.size() < 2) return false;

			// 投影到主轴，获取 t 值
			int axis = 0; // 选择 D 最大分量作为主轴以保证精度
			if (std::abs(D.y()) > std::abs(D.x())) axis = 1;
			if (std::abs(D.z()) > std::abs(D.x()) && std::abs(D.z()) > std::abs(D.y())) axis = 2;

			double v0 = (axis == 0 ? pts[0].x() : (axis == 1 ? pts[0].y() : pts[0].z()));
			double v1 = (axis == 0 ? pts[1].x() : (axis == 1 ? pts[1].y() : pts[1].z()));

			if (v0 > v1) std::swap(v0, v1);
			t_min = v0; t_max = v1;
			return true;
			};

		double t1_min, t1_max, t2_min, t2_max;
		if (!get_inter_interval(t1, dv0, dv1, dv2, t1_min, t1_max)) return false;
		if (!get_inter_interval(t2, du0, du1, du2, t2_min, t2_max)) return false;

		// 检查区间重叠
		double t_start = std::max(t1_min, t2_min);
		double t_end = std::min(t1_max, t2_max);

		if (t_start > t_end) return false; // 不重叠

		// 反算回 3D 坐标
		// 我们需要 Line 上的一个基点，这里用 n1, n2, d1, d2 求解
		// 简化：直接插值回 3D
		// 我们知道 t 是 D 向量分量上的投影。
		// 为了简便，我们直接判定：相交线段的端点是有效的
		// 更准确的方法是求解 P_start = Origin + t_start * D (需确定 Origin)

		// 替代方案：直接计算四个交点，取中间两个
		// 这种情况下，直接用区间比例回溯比较麻烦。
		// 我们用一种更简单的方法：既然已经确认相交，
		// 我们直接找 (T1 截 Plane2 的线段) 和 (T2 截 Plane1 的线段) 在 3D 空间是否有交集（共线）

		// 实际上，t_start 和 t_end 对应的就是在 D 上的投影值。
		// 我们只需找到 Line 上的一个点 P0，那么 Intersection = P0 + t * D ? 不完全是，取决于投影轴。
		// 让我们用最笨但有效的方法：3D 坐标插值

		// 重新计算线段端点 (T1 in Plane2)
		std::vector<SPAposition> seg1;
		if (dv0 * dv1 <= EPSILON) seg1.push_back(t1.p0 + (t1.p1 - t1.p0) * (-dv0 / (dv1 - dv0 + 1e-15)));
		if (dv1 * dv2 <= EPSILON) seg1.push_back(t1.p1 + (t1.p2 - t1.p1) * (-dv1 / (dv2 - dv1 + 1e-15)));
		if (dv2 * dv0 <= EPSILON) seg1.push_back(t1.p2 + (t1.p0 - t1.p2) * (-dv2 / (dv0 - dv2 + 1e-15)));

		// 重新计算线段端点 (T2 in Plane1)
		std::vector<SPAposition> seg2;
		if (du0 * du1 <= EPSILON) seg2.push_back(t2.p0 + (t2.p1 - t2.p0) * (-du0 / (du1 - du0 + 1e-15)));
		if (du1 * du2 <= EPSILON) seg2.push_back(t2.p1 + (t2.p2 - t2.p1) * (-du1 / (du2 - du1 + 1e-15)));
		if (du2 * du0 <= EPSILON) seg2.push_back(t2.p2 + (t2.p0 - t2.p2) * (-du2 / (du0 - du2 + 1e-15)));

		if (seg1.size() < 2 || seg2.size() < 2) return false;

		// 现在我们在同一条直线上有线段 S1(a,b) 和 S2(c,d)
		// 求它们的交集
		// 投影到 D 上排序
		struct Node { double t; SPAposition p; int id; }; // id: 0=S1, 1=S2
		std::vector<Node> nodes;

		auto add_nodes = [&](const std::vector<SPAposition>& seg, int id) {
			double val0 = SPAposition::dot(seg[0], D);
			double val1 = SPAposition::dot(seg[1], D);
			if (val0 < val1) { nodes.push_back({ val0, seg[0], id }); nodes.push_back({ val1, seg[1], id }); }
			else { nodes.push_back({ val1, seg[1], id }); nodes.push_back({ val0, seg[0], id }); }
			};

		add_nodes(seg1, 0);
		add_nodes(seg2, 1);

		double max_min = std::max(nodes[0].t, nodes[2].t); // S1_start, S2_start
		double min_max = std::min(nodes[1].t, nodes[3].t); // S1_end, S2_end

		if (max_min > min_max + EPSILON) return false; // 无重叠

		// 重叠部分的端点
		// 找到对应 t 值的点有点麻烦，简单起见，取几何上位于中间的两个点
		// 排序所有4个点在D上的投影
		std::vector<std::pair<double, SPAposition>> sort_pts;
		for (auto& p : seg1) sort_pts.push_back({ SPAposition::dot(p, D), p });
		for (auto& p : seg2) sort_pts.push_back({ SPAposition::dot(p, D), p });
		std::sort(sort_pts.begin(), sort_pts.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

		out1 = sort_pts[1].second; // 中间两个点即为交集
		out2 = sort_pts[2].second;

		return true;
	}
}

// ==========================================
// 3. 数据结构部分 (保持原样)
// ==========================================
struct TriPtr {
	int surface_id; int type; int u_idx; int v_idx;
	bool operator<(const TriPtr& o) const { return std::tie(surface_id, type, u_idx, v_idx) < std::tie(o.surface_id, o.type, o.u_idx, o.v_idx); }
};

class Surface {
public:
	std::vector<SPAposition> points;
	int u_count, v_count, id;
	Surface(int u, int v, int pid) : u_count(u), v_count(v), id(pid) { points.resize(u * v); }
	void setPoint(int u, int v, double x, double y, double z) { if (u >= 0 && u < u_count && v >= 0 && v < v_count) points[v * u_count + u] = SPAposition(x, y, z); }
	SPAposition getPoint(int u, int v) const { return (u >= 0 && u < u_count && v >= 0 && v < v_count) ? points[v * u_count + u] : SPAposition(0, 0, 0); }

	Triangle getTriangleCoords(const TriPtr& tri) const {
		int u = tri.u_idx, v = tri.v_idx;
		if (tri.type == 1) return { getPoint(u, v), getPoint(u + 1, v), getPoint(u, v + 1) };
		else return { getPoint(u + 1, v), getPoint(u + 1, v + 1), getPoint(u, v + 1) };
	}
	SPAposition getCentroid(const TriPtr& tri) const {
		Triangle t = getTriangleCoords(tri);
		return (t.p0 + t.p1 + t.p2) / 3.0;
	}
	bool isValid(const TriPtr& tri) const { return (tri.u_idx >= 0 && tri.u_idx < u_count - 1 && tri.v_idx >= 0 && tri.v_idx < v_count - 1); }
};

// ==========================================
// 4. Solver (逻辑修改)
// ==========================================
class Solver {
private:
	Surface& surfA; Surface& surfB; std::ofstream& logFile;
	std::set<TriPtr> global_visited_A, global_visited_B;

public:
	Solver(Surface& a, Surface& b, std::ofstream& fs) : surfA(a), surfB(b), logFile(fs) {}

	double calcCentroidDistSq(const TriPtr& ta, const TriPtr& tb) {
		return (surfA.getCentroid(ta) - surfB.getCentroid(tb)).magSq();
	}

	std::vector<TriPtr> getNextNeighbors(const TriPtr& current, const Surface& surf) {
		std::vector<TriPtr> neighbors;
		int u = current.u_idx; int v = current.v_idx;
		if (current.type == 1) {
			TriPtr n = { current.surface_id, 2, u, v };
			if (surf.isValid(n)) neighbors.push_back(n);
		}
		else {
			TriPtr n1 = { current.surface_id, 1, u + 1, v }; if (surf.isValid(n1)) neighbors.push_back(n1);
			TriPtr n2 = { current.surface_id, 1, u, v + 1 }; if (surf.isValid(n2)) neighbors.push_back(n2);
		}
		return neighbors;
	}

	void run() {
		logFile << "SECTION: PROCESS" << std::endl;
		for (int va = 0; va < surfA.v_count - 1; ++va) {
			for (int ua = 0; ua < surfA.u_count - 1; ++ua) {
				TriPtr startA = { 0, 1, ua, va };
				if (!surfA.isValid(startA) || global_visited_A.count(startA)) continue;

				for (int vb = 0; vb < surfB.v_count - 1; ++vb) {
					for (int ub = 0; ub < surfB.u_count - 1; ++ub) {
						TriPtr startB = { 1, 1, ub, vb };
						if (!surfB.isValid(startB) || global_visited_B.count(startB)) continue;
						if (global_visited_A.count(startA)) break;

						processPath(startA, startB);
					}
				}
			}
		}
	}

	void processPath(TriPtr currA, TriPtr currB) {
		logFile << "PATH_START" << std::endl;
		bool alive = true;
		while (alive) {
			global_visited_A.insert(currA);
			global_visited_B.insert(currB);

			double distSq = calcCentroidDistSq(currA, currB);
			logFile << "STEP A " << currA.type << " " << currA.u_idx << " " << currA.v_idx << " "
				<< "B " << currB.type << " " << currB.u_idx << " " << currB.v_idx << " "
				<< "Dist " << std::sqrt(distSq) << std::endl;

			auto nextAs = getNextNeighbors(currA, surfA);
			auto nextBs = getNextNeighbors(currB, surfB);

			double bestDist = distSq;
			TriPtr nextA = currA; TriPtr nextB = currB;
			bool moved = false;

			for (const auto& na : nextAs) {
				double d = calcCentroidDistSq(na, currB);
				if (d < bestDist) { bestDist = d; nextA = na; nextB = currB; moved = true; }
			}
			for (const auto& nb : nextBs) {
				double d = calcCentroidDistSq(currA, nb);
				if (d < bestDist) { bestDist = d; nextA = currA; nextB = nb; moved = true; }
			}

			if (moved) {
				currA = nextA; currB = nextB;
			}
			else {
				outputResult(currA, currB);
				logFile << "END DEAD_POINTER" << std::endl;
				alive = false;
			}
		}
	}

	// --- 最终输出逻辑：区分相交与最近点 ---
	void outputResult(const TriPtr& ta, const TriPtr& tb) {
		Triangle triA = surfA.getTriangleCoords(ta);
		Triangle triB = surfB.getTriangleCoords(tb);
		SPAposition i1, i2;

		// 1. 优先检测是否相交
		if (Geometry::intersectTriTri(triA, triB, i1, i2)) {
			// 输出交线
			logFile << "RESULT_INTERSECTION "
				<< i1.x() << " " << i1.y() << " " << i1.z() << " "
				<< i2.x() << " " << i2.y() << " " << i2.z() << " "
				<< "0.0" << std::endl;
		}
		else {
			// 2. 如果不相交，计算最近点对
			SPAposition bestA, bestB;
			double minDistSq = 1e30;

			// A 顶点 -> B 面
			SPAposition vertsA[3] = { triA.p0, triA.p1, triA.p2 };
			for (auto& p : vertsA) {
				SPAposition c = Geometry::closestPointOnTriangle(p, triB);
				double d = (p - c).magSq();
				if (d < minDistSq) { minDistSq = d; bestA = p; bestB = c; }
			}
			// B 顶点 -> A 面
			SPAposition vertsB[3] = { triB.p0, triB.p1, triB.p2 };
			for (auto& p : vertsB) {
				SPAposition c = Geometry::closestPointOnTriangle(p, triA);
				double d = (p - c).magSq();
				if (d < minDistSq) { minDistSq = d; bestB = p; bestA = c; }
			}

			logFile << "RESULT_NEAREST "
				<< bestA.x() << " " << bestA.y() << " " << bestA.z() << " "
				<< bestB.x() << " " << bestB.y() << " " << bestB.z() << " "
				<< std::sqrt(minDistSq) << std::endl;
		}
	}
};

// ==========================================
// 5. Main
// ==========================================
void initData(Surface& A, Surface& B) {
	// Surface A
	A.setPoint(0, 0, -1.0, -1.0, 0.0); A.setPoint(1, 0, -0.5, -1.0, 0.1); A.setPoint(2, 0, 0.0, -1.0, 0.0); A.setPoint(3, 0, 0.5, -1.0, 0.1); A.setPoint(4, 0, 1.0, -1.0, 0.0);
	A.setPoint(0, 1, -1.0, -0.5, 0.1); A.setPoint(1, 1, -0.5, -0.5, 0.3); A.setPoint(2, 1, 0.0, -0.5, 0.2); A.setPoint(3, 1, 0.5, -0.5, 0.3); A.setPoint(4, 1, 1.0, -0.5, 0.1);
	A.setPoint(0, 2, -1.0, 0.0, 0.0);  A.setPoint(1, 2, -0.5, 0.0, 0.2);  A.setPoint(2, 2, 0.0, 0.0, 0.5);  A.setPoint(3, 2, 0.5, 0.0, 0.2);  A.setPoint(4, 2, 1.0, 0.0, 0.0);
	A.setPoint(0, 3, -1.0, 0.5, 0.1);  A.setPoint(1, 3, -0.5, 0.5, 0.3);  A.setPoint(2, 3, 0.0, 0.5, 0.2);  A.setPoint(3, 3, 0.5, 0.5, 0.3);  A.setPoint(4, 3, 1.0, 0.5, 0.1);
	A.setPoint(0, 4, -1.0, 1.0, 0.0);  A.setPoint(1, 4, -0.5, 1.0, 0.1);  A.setPoint(2, 4, 0.0, 1.0, 0.0);  A.setPoint(3, 4, 0.5, 1.0, 0.1);  A.setPoint(4, 4, 1.0, 1.0, 0.0);

	// Surface B (Constructed to intersect A)
	B.setPoint(0, 0, -1.0, -1.0, 0.5); B.setPoint(1, 0, -0.5, -1.0, 0.4); B.setPoint(2, 0, 0.0, -1.0, 0.5); B.setPoint(3, 0, 0.5, -1.0, 0.4); B.setPoint(4, 0, 1.0, -1.0, 0.5);
	B.setPoint(0, 1, -1.0, -0.5, 0.4); B.setPoint(1, 1, -0.5, -0.5, 0.1); B.setPoint(2, 1, 0.0, -0.5, -0.2); B.setPoint(3, 1, 0.5, -0.5, 0.1); B.setPoint(4, 1, 1.0, -0.5, 0.4);
	B.setPoint(0, 2, -1.0, 0.0, 0.5);  B.setPoint(1, 2, -0.5, 0.0, -0.2); B.setPoint(2, 2, 0.0, 0.0, -1.0); B.setPoint(3, 2, 0.5, 0.0, -0.2); B.setPoint(4, 2, 1.0, 0.0, 0.5);
	B.setPoint(0, 3, -1.0, 0.5, 0.4);  B.setPoint(1, 3, -0.5, 0.5, 0.1);  B.setPoint(2, 3, 0.0, 0.5, -0.2); B.setPoint(3, 3, 0.5, 0.5, 0.1); B.setPoint(4, 3, 1.0, 0.5, 0.4);
	B.setPoint(0, 4, -1.0, 1.0, 0.0);  B.setPoint(1, 4, -0.5, 1.0, 0.4);  B.setPoint(2, 4, 0.0, 1.0, 0.5);  B.setPoint(3, 4, 0.5, 1.0, 0.4);  B.setPoint(4, 4, 1.0, 1.0, 0.5);
}

void exportGeometry(const Surface& surf, const std::string& name, std::ofstream& fs) {
	fs << "SECTION: GEOMETRY " << name << std::endl;
	fs << "DIMS " << surf.u_count << " " << surf.v_count << std::endl;
	for (int v = 0; v < surf.v_count; ++v) {
		for (int u = 0; u < surf.u_count; ++u) {
			SPAposition p = surf.getPoint(u, v);
			fs << p.x() << " " << p.y() << " " << p.z() << std::endl;
		}
	}
}

int main() {
	Surface surfA(5, 5, 0);
	Surface surfB(5, 5, 1);
	initData(surfA, surfB);

	std::ofstream fs("debug_output.txt");
	if (!fs.is_open()) return 1;

	exportGeometry(surfA, "A", fs);
	exportGeometry(surfB, "B", fs);

	Solver solver(surfA, surfB, fs);
	solver.run();

	fs.close();
	std::cout << "Complete. Intersection lines and nearest pairs calculated." << std::endl;
	return 0;
}