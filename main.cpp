#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <tuple>
#include <vector>

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
};

struct Triangle {
    SPAposition p0, p1, p2;
};

namespace Geometry {

SPAposition closestPointOnSegment(const SPAposition& p, const SPAposition& a, const SPAposition& b) {
    SPAposition ab = b - a;
    double denom = std::max(ab.magSq(), EPSILON);
    double t = SPAposition::dot(p - a, ab) / denom;
    t = std::max(0.0, std::min(1.0, t));
    return a + ab * t;
}

SPAposition closestPointOnTriangle(const SPAposition& p, const Triangle& tri) {
    SPAposition ab = tri.p1 - tri.p0;
    SPAposition ac = tri.p2 - tri.p0;
    SPAposition ap = p - tri.p0;

    SPAposition n = SPAposition::cross(ab, ac);
    double nSq = std::max(n.magSq(), EPSILON);
    SPAposition p_proj = p - n * (SPAposition::dot(ap, n) / nSq);

    SPAposition bc = tri.p2 - tri.p1;
    SPAposition ca = tri.p0 - tri.p2;

    if (SPAposition::dot(SPAposition::cross(ab, n), p_proj - tri.p0) > 0) goto check_edges;
    if (SPAposition::dot(SPAposition::cross(bc, n), p_proj - tri.p1) > 0) goto check_edges;
    if (SPAposition::dot(SPAposition::cross(ca, n), p_proj - tri.p2) > 0) goto check_edges;

    return p_proj;

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

bool intersectTriTri(const Triangle& t1, const Triangle& t2, SPAposition& out1, SPAposition& out2) {
    SPAposition n1 = SPAposition::cross(t1.p1 - t1.p0, t1.p2 - t1.p0);
    SPAposition n2 = SPAposition::cross(t2.p1 - t2.p0, t2.p2 - t2.p0);
    SPAposition d = SPAposition::cross(n1, n2);

    if (d.magSq() < EPSILON) return false;

    auto signedDist = [](const SPAposition& n, const SPAposition& p0, const SPAposition& p) {
        return SPAposition::dot(n, p - p0);
    };

    double du0 = signedDist(n1, t1.p0, t2.p0), du1 = signedDist(n1, t1.p0, t2.p1), du2 = signedDist(n1, t1.p0, t2.p2);
    double dv0 = signedDist(n2, t2.p0, t1.p0), dv1 = signedDist(n2, t2.p0, t1.p1), dv2 = signedDist(n2, t2.p0, t1.p2);

    auto sameSide = [](double a, double b, double c) {
        if (std::abs(a) < EPSILON || std::abs(b) < EPSILON || std::abs(c) < EPSILON) return false;
        return (a > 0 && b > 0 && c > 0) || (a < 0 && b < 0 && c < 0);
    };

    if (sameSide(du0, du1, du2) || sameSide(dv0, dv1, dv2)) return false;

    auto collectSegment = [](const Triangle& tri, double d0, double d1, double d2) {
        std::vector<SPAposition> pts;
        auto addIfCross = [&](const SPAposition& a, const SPAposition& b, double da, double db) {
            if (da * db > EPSILON) return;
            double denom = db - da;
            if (std::abs(denom) < EPSILON) return;
            double t = -da / denom;
            pts.push_back(a + (b - a) * t);
        };
        addIfCross(tri.p0, tri.p1, d0, d1);
        addIfCross(tri.p1, tri.p2, d1, d2);
        addIfCross(tri.p2, tri.p0, d2, d0);
        return pts;
    };

    auto seg1 = collectSegment(t1, dv0, dv1, dv2);
    auto seg2 = collectSegment(t2, du0, du1, du2);
    if (seg1.size() < 2 || seg2.size() < 2) return false;

    std::vector<std::pair<double, SPAposition>> pts;
    for (const auto& p : seg1) pts.push_back({SPAposition::dot(p, d), p});
    for (const auto& p : seg2) pts.push_back({SPAposition::dot(p, d), p});
    std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    out1 = pts[1].second;
    out2 = pts[2].second;
    return true;
}

} // namespace Geometry

struct TriPtr {
    int surface_id;
    int type;
    int u_idx;
    int v_idx;

    bool operator<(const TriPtr& o) const {
        return std::tie(surface_id, type, u_idx, v_idx) < std::tie(o.surface_id, o.type, o.u_idx, o.v_idx);
    }
};

class Surface {
public:
    std::vector<SPAposition> points;
    int u_count, v_count, id;

    Surface(int u, int v, int pid) : u_count(u), v_count(v), id(pid), points(static_cast<size_t>(u * v)) {}

    void setPoint(int u, int v, double x, double y, double z) {
        if (u >= 0 && u < u_count && v >= 0 && v < v_count) points[v * u_count + u] = SPAposition(x, y, z);
    }

    SPAposition getPoint(int u, int v) const {
        if (u < 0 || u >= u_count || v < 0 || v >= v_count) return SPAposition();
        return points[v * u_count + u];
    }

    bool isValid(const TriPtr& tri) const {
        return tri.u_idx >= 0 && tri.u_idx < u_count - 1 && tri.v_idx >= 0 && tri.v_idx < v_count - 1 && (tri.type == 1 || tri.type == 2);
    }

    Triangle getTriangleCoords(const TriPtr& tri) const {
        int u = tri.u_idx, v = tri.v_idx;
        if (tri.type == 1) return {getPoint(u, v), getPoint(u + 1, v), getPoint(u, v + 1)};
        return {getPoint(u + 1, v), getPoint(u + 1, v + 1), getPoint(u, v + 1)};
    }

    SPAposition getCentroid(const TriPtr& tri) const {
        Triangle t = getTriangleCoords(tri);
        return (t.p0 + t.p1 + t.p2) / 3.0;
    }

    std::vector<TriPtr> allTriangles() const {
        std::vector<TriPtr> tris;
        tris.reserve(static_cast<size_t>((u_count - 1) * (v_count - 1) * 2));
        for (int v = 0; v < v_count - 1; ++v) {
            for (int u = 0; u < u_count - 1; ++u) {
                tris.push_back({id, 1, u, v});
                tris.push_back({id, 2, u, v});
            }
        }
        return tris;
    }
};

struct RunStats {
    int path_count = 0;
    int step_count = 0;
};

class Solver {
private:
    Surface& surfA;
    Surface& surfB;
    std::ofstream& logFile;
    std::set<TriPtr> visitedA, visitedB;

public:
    Solver(Surface& a, Surface& b, std::ofstream& fs) : surfA(a), surfB(b), logFile(fs) {}

    double centroidDistSq(const TriPtr& ta, const TriPtr& tb) {
        return (surfA.getCentroid(ta) - surfB.getCentroid(tb)).magSq();
    }

    std::vector<TriPtr> nextNeighbors(const TriPtr& tri, const Surface& surf) {
        std::vector<TriPtr> n;
        int u = tri.u_idx, v = tri.v_idx;

        TriPtr sameCell = {tri.surface_id, tri.type == 1 ? 2 : 1, u, v};
        TriPtr right = {tri.surface_id, tri.type, u + 1, v};
        TriPtr up = {tri.surface_id, tri.type, u, v + 1};

        if (surf.isValid(sameCell)) n.push_back(sameCell);
        if (surf.isValid(right)) n.push_back(right);
        if (surf.isValid(up)) n.push_back(up);

        return n;
    }

    void outputResult(const TriPtr& ta, const TriPtr& tb) {
        Triangle triA = surfA.getTriangleCoords(ta);
        Triangle triB = surfB.getTriangleCoords(tb);

        SPAposition i1, i2;
        if (Geometry::intersectTriTri(triA, triB, i1, i2)) {
            logFile << "RESULT_INTERSECTION " << i1.x() << " " << i1.y() << " " << i1.z() << " " << i2.x() << " " << i2.y() << " " << i2.z() << " 0.0\n";
            return;
        }

        SPAposition bestA, bestB;
        double minDistSq = std::numeric_limits<double>::max();
        SPAposition vertsA[3] = {triA.p0, triA.p1, triA.p2};
        SPAposition vertsB[3] = {triB.p0, triB.p1, triB.p2};

        for (const auto& p : vertsA) {
            SPAposition c = Geometry::closestPointOnTriangle(p, triB);
            double d = (p - c).magSq();
            if (d < minDistSq) {
                minDistSq = d;
                bestA = p;
                bestB = c;
            }
        }
        for (const auto& p : vertsB) {
            SPAposition c = Geometry::closestPointOnTriangle(p, triA);
            double d = (p - c).magSq();
            if (d < minDistSq) {
                minDistSq = d;
                bestA = c;
                bestB = p;
            }
        }

        logFile << "RESULT_NEAREST " << bestA.x() << " " << bestA.y() << " " << bestA.z() << " " << bestB.x() << " " << bestB.y() << " " << bestB.z() << " " << std::sqrt(minDistSq) << "\n";
    }

    void processPath(TriPtr a, TriPtr b, RunStats& stats) {
        logFile << "PATH_START\n";
        stats.path_count += 1;

        while (true) {
            visitedA.insert(a);
            visitedB.insert(b);
            stats.step_count += 1;

            double cur = centroidDistSq(a, b);
            logFile << "STEP A " << a.type << " " << a.u_idx << " " << a.v_idx << " B " << b.type << " " << b.u_idx << " " << b.v_idx << " Dist " << std::sqrt(cur) << "\n";

            TriPtr bestA = a, bestB = b;
            double best = cur;

            for (const auto& na : nextNeighbors(a, surfA)) {
                double d = centroidDistSq(na, b);
                if (d + EPSILON < best) {
                    best = d;
                    bestA = na;
                    bestB = b;
                }
            }
            for (const auto& nb : nextNeighbors(b, surfB)) {
                double d = centroidDistSq(a, nb);
                if (d + EPSILON < best) {
                    best = d;
                    bestA = a;
                    bestB = nb;
                }
            }

            if (bestA.u_idx == a.u_idx && bestA.v_idx == a.v_idx && bestA.type == a.type &&
                bestB.u_idx == b.u_idx && bestB.v_idx == b.v_idx && bestB.type == b.type) {
                outputResult(a, b);
                logFile << "END DEAD_POINTER\n";
                break;
            }

            a = bestA;
            b = bestB;
        }
    }

    RunStats run() {
        RunStats stats;
        logFile << "SECTION: PROCESS\n";

        for (int va = 0; va < surfA.v_count - 1; ++va) {
            for (int ua = 0; ua < surfA.u_count - 1; ++ua) {
                for (int ta : {1, 2}) {
                    TriPtr startA = {0, ta, ua, va};
                    if (!surfA.isValid(startA) || visitedA.count(startA)) continue;

                    for (int vb = 0; vb < surfB.v_count - 1; ++vb) {
                        for (int ub = 0; ub < surfB.u_count - 1; ++ub) {
                            for (int tb : {1, 2}) {
                                TriPtr startB = {1, tb, ub, vb};
                                if (!surfB.isValid(startB) || visitedB.count(startB)) continue;
                                if (visitedA.count(startA)) break;
                                processPath(startA, startB, stats);
                            }
                        }
                    }
                }
            }
        }

        return stats;
    }
};

void exportGeometry(const Surface& surf, const std::string& name, std::ofstream& fs) {
    fs << "SECTION: GEOMETRY " << name << "\n";
    fs << "DIMS " << surf.u_count << " " << surf.v_count << "\n";
    for (int v = 0; v < surf.v_count; ++v) {
        for (int u = 0; u < surf.u_count; ++u) {
            SPAposition p = surf.getPoint(u, v);
            fs << p.x() << " " << p.y() << " " << p.z() << "\n";
        }
    }
}

void initCaseData(Surface& A, Surface& B, int case_id) {
    for (int v = 0; v < A.v_count; ++v) {
        for (int u = 0; u < A.u_count; ++u) {
            double x = -1.2 + 2.4 * static_cast<double>(u) / (A.u_count - 1);
            double y = -1.2 + 2.4 * static_cast<double>(v) / (A.v_count - 1);

            double zA = 0.25 * std::sin((case_id + 1) * x) + 0.20 * std::cos((case_id + 2) * y);
            double zB = 0.22 * std::cos((case_id + 2) * x) - 0.25 * std::sin((case_id + 1) * y);

            if (case_id % 2 == 0) {
                zB += 0.08 * x - 0.05;
            } else {
                zA += 0.05 * y;
                zB -= 0.03 * x;
            }

            if (case_id == 4) zB -= 0.18;
            if (case_id == 5) zA += 0.15;

            A.setPoint(u, v, x, y, zA);
            B.setPoint(u, v, x, y, zB);
        }
    }
}

struct BaselineStats {
    double best_dist = 0.0;
    long long pair_count = 0;
};

BaselineStats bruteForceCentroid(const Surface& A, const Surface& B) {
    auto trisA = A.allTriangles();
    auto trisB = B.allTriangles();

    BaselineStats s;
    s.pair_count = static_cast<long long>(trisA.size()) * static_cast<long long>(trisB.size());

    double best = std::numeric_limits<double>::max();
    for (const auto& ta : trisA) {
        SPAposition ca = A.getCentroid(ta);
        for (const auto& tb : trisB) {
            SPAposition cb = B.getCentroid(tb);
            best = std::min(best, (ca - cb).magSq());
        }
    }

    s.best_dist = std::sqrt(best);
    return s;
}

int main() {
    std::filesystem::create_directories("outputs");

    const int case_count = 6;
    const int grid_u = 8;
    const int grid_v = 8;

    for (int case_id = 0; case_id < case_count; ++case_id) {
        Surface surfA(grid_u, grid_v, 0);
        Surface surfB(grid_u, grid_v, 1);
        initCaseData(surfA, surfB, case_id);

        std::string output_file = "outputs/case_" + std::to_string(case_id + 1) + "_debug_output.txt";
        std::ofstream fs(output_file);
        if (!fs.is_open()) {
            std::cerr << "Failed to open " << output_file << "\n";
            return 1;
        }

        fs << "SECTION: META\n";
        fs << "CASE_ID " << (case_id + 1) << "\n";

        exportGeometry(surfA, "A", fs);
        exportGeometry(surfB, "B", fs);

        auto greedy_start = std::chrono::steady_clock::now();
        Solver solver(surfA, surfB, fs);
        RunStats run_stats = solver.run();
        auto greedy_end = std::chrono::steady_clock::now();

        auto brute_start = std::chrono::steady_clock::now();
        BaselineStats brute_stats = bruteForceCentroid(surfA, surfB);
        auto brute_end = std::chrono::steady_clock::now();

        double greedy_ms = std::chrono::duration<double, std::milli>(greedy_end - greedy_start).count();
        double brute_ms = std::chrono::duration<double, std::milli>(brute_end - brute_start).count();

        fs << "SECTION: METRICS\n";
        fs << "METRIC GREEDY_MS " << greedy_ms << "\n";
        fs << "METRIC BRUTE_FORCE_MS " << brute_ms << "\n";
        fs << "METRIC GREEDY_PATHS " << run_stats.path_count << "\n";
        fs << "METRIC GREEDY_STEPS " << run_stats.step_count << "\n";
        fs << "METRIC BRUTE_FORCE_PAIRS " << brute_stats.pair_count << "\n";
        fs << "METRIC BRUTE_FORCE_BEST_DIST " << brute_stats.best_dist << "\n";

        std::cout << "Generated " << output_file << "\n";
    }

    std::cout << "Complete. 6 cases exported to outputs/." << std::endl;
    return 0;
}
