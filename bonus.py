import glob
import os
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def parse_case_file(filename):
    with open(filename, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    data_A, data_B = [], []
    paths = []
    metrics = {}
    case_id = None

    current_section = None
    dims = (0, 0)
    temp_points = []
    current_path = []

    for line in lines:
        if line.startswith("SECTION: META"):
            current_section = "META"
        elif line.startswith("CASE_ID"):
            case_id = int(line.split()[1])
        elif line.startswith("SECTION: GEOMETRY"):
            current_section = "A" if line.endswith("A") else "B"
            temp_points = []
        elif line.startswith("DIMS"):
            dims = tuple(map(int, line.split()[1:]))
        elif line.startswith("SECTION: PROCESS"):
            current_section = "PROCESS"
        elif line.startswith("SECTION: METRICS"):
            current_section = "METRICS"
        elif line.startswith("PATH_START"):
            current_path = []
        elif line.startswith("END"):
            if current_path:
                paths.append(current_path)
            current_path = []
        elif line.startswith("METRIC"):
            _, key, value = line.split(maxsplit=2)
            try:
                metrics[key] = float(value)
            except ValueError:
                metrics[key] = value
        elif line.startswith("RESULT_INTERSECTION") or line.startswith("RESULT_NEAREST"):
            parts = line.split()
            current_path.append(
                {
                    "type": "result",
                    "res_type": "intersection" if "INTERSECTION" in line else "nearest",
                    "pA": [float(parts[1]), float(parts[2]), float(parts[3])],
                    "pB": [float(parts[4]), float(parts[5]), float(parts[6])],
                    "dist": float(parts[7]),
                }
            )
        elif line.startswith("STEP"):
            parts = line.split()
            current_path.append(
                {
                    "type": "step",
                    "A": {"type": int(parts[2]), "u": int(parts[3]), "v": int(parts[4])},
                    "B": {"type": int(parts[6]), "u": int(parts[7]), "v": int(parts[8])},
                    "dist": float(parts[10]),
                }
            )
        else:
            if current_section in ["A", "B"]:
                coords = list(map(float, line.split()))
                temp_points.append(coords)
                if dims[0] > 0 and len(temp_points) == dims[0] * dims[1]:
                    grid = np.array(temp_points).reshape(dims[1], dims[0], 3)
                    if current_section == "A":
                        data_A = grid
                    else:
                        data_B = grid

    return {
        "case_id": case_id,
        "file": filename,
        "A": data_A,
        "B": data_B,
        "paths": paths,
        "metrics": metrics,
    }


def create_mesh_wireframe(grid, color, name, opacity=0.35):
    x, y, z = [], [], []
    rows, cols, _ = grid.shape
    for r in range(rows):
        for c in range(cols):
            x.append(grid[r, c, 0])
            y.append(grid[r, c, 1])
            z.append(grid[r, c, 2])
        x.append(None)
        y.append(None)
        z.append(None)
    for c in range(cols):
        for r in range(rows):
            x.append(grid[r, c, 0])
            y.append(grid[r, c, 1])
            z.append(grid[r, c, 2])
        x.append(None)
        y.append(None)
        z.append(None)
    return go.Scatter3d(x=x, y=y, z=z, mode="lines", line=dict(color=color, width=2), opacity=opacity, name=name)


def collect_result_lines(paths):
    inter_x, inter_y, inter_z = [], [], []
    near_x, near_y, near_z = [], [], []
    for path in paths:
        res = next((item for item in path if item["type"] == "result"), None)
        if not res:
            continue
        if res["res_type"] == "intersection":
            inter_x.extend([res["pA"][0], res["pB"][0], None])
            inter_y.extend([res["pA"][1], res["pB"][1], None])
            inter_z.extend([res["pA"][2], res["pB"][2], None])
        else:
            near_x.extend([res["pA"][0], res["pB"][0], None])
            near_y.extend([res["pA"][1], res["pB"][1], None])
            near_z.extend([res["pA"][2], res["pB"][2], None])

    trace_inter = go.Scatter3d(x=inter_x, y=inter_y, z=inter_z, mode="lines", line=dict(color="red", width=8), name="Intersection")
    trace_near = go.Scatter3d(x=near_x, y=near_y, z=near_z, mode="lines", line=dict(color="gold", width=5), name="Nearest")
    return trace_inter, trace_near


def create_case_summary(case_data, output_dir):
    wire_a = create_mesh_wireframe(case_data["A"], "royalblue", "A")
    wire_b = create_mesh_wireframe(case_data["B"], "orangered", "B")
    trace_inter, trace_near = collect_result_lines(case_data["paths"])

    fig = go.Figure(data=[wire_a, wire_b, trace_inter, trace_near])
    fig.update_layout(
        title=f"Case {case_data['case_id']} | Search Results",
        scene=dict(aspectmode="data"),
        legend=dict(x=0.02, y=0.98),
    )
    output_file = output_dir / f"case_{case_data['case_id']}_summary.html"
    fig.write_html(output_file)


def create_timing_dashboard(case_items, output_dir):
    case_labels = [f"Case {c['case_id']}" for c in case_items]
    greedy_ms = [c["metrics"].get("GREEDY_MS", 0.0) for c in case_items]
    brute_ms = [c["metrics"].get("BRUTE_FORCE_MS", 0.0) for c in case_items]
    steps = [c["metrics"].get("GREEDY_STEPS", 0.0) for c in case_items]
    pairs = [c["metrics"].get("BRUTE_FORCE_PAIRS", 0.0) for c in case_items]

    speedup = [(b / g) if g > 1e-9 else 0 for b, g in zip(brute_ms, greedy_ms)]

    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=("Time Cost Comparison", "Speedup (Brute / Greedy)", "Greedy Steps", "Brute-force Pair Count"),
    )

    fig.add_trace(go.Bar(x=case_labels, y=greedy_ms, name="My Search (Greedy)", marker_color="teal"), row=1, col=1)
    fig.add_trace(go.Bar(x=case_labels, y=brute_ms, name="Brute-force Pairwise", marker_color="tomato"), row=1, col=1)

    fig.add_trace(go.Scatter(x=case_labels, y=speedup, mode="lines+markers", name="Speedup", line=dict(color="purple", width=3)), row=1, col=2)
    fig.add_trace(go.Bar(x=case_labels, y=steps, name="Greedy Steps", marker_color="deepskyblue"), row=2, col=1)
    fig.add_trace(go.Bar(x=case_labels, y=pairs, name="Brute Pairs", marker_color="gray"), row=2, col=2)

    fig.update_layout(title="Greedy vs Brute-force Benchmark Dashboard", barmode="group", height=900)
    fig.update_yaxes(title_text="Milliseconds", row=1, col=1)
    fig.update_yaxes(title_text="x", row=1, col=2)

    output_file = output_dir / "timing_dashboard.html"
    fig.write_html(output_file)


def main():
    files = sorted(glob.glob("outputs/case_*_debug_output.txt"))
    if not files:
        print("No case files found under outputs/. Please run C++ first.")
        return

    case_items = [parse_case_file(f) for f in files]
    case_items = sorted(case_items, key=lambda x: x["case_id"])

    output_dir = Path("outputs")
    for item in case_items:
        create_case_summary(item, output_dir)

    create_timing_dashboard(case_items, output_dir)

    print(f"Loaded {len(case_items)} cases.")
    print("Generated case summary HTML files and outputs/timing_dashboard.html")


if __name__ == "__main__":
    main()
