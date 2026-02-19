import plotly.graph_objects as go
import numpy as np
import os


# ==========================================
# 1. 核心解析器
# ==========================================
def parse_file(filename):
    # 优先尝试指定路径，如果不存在则尝试当前目录
    if not os.path.exists(filename):
        local_filename = os.path.basename(filename)
        if os.path.exists(local_filename):
            filename = local_filename
        else:
            print(f"Error: File {filename} not found.")
            return None, None, []

    with open(filename, 'r') as f:
        lines = f.readlines()

    data_A = []
    data_B = []
    paths = []

    current_section = None
    temp_points = []
    dims = (0, 0)
    current_path = []

    for line in lines:
        line = line.strip()
        if not line: continue

        if line.startswith("SECTION: GEOMETRY"):
            if "A" in line:
                current_section = "A"
            else:
                current_section = "B"
            temp_points = []
        elif line.startswith("DIMS"):
            dims = tuple(map(int, line.split()[1:]))
        elif line.startswith("SECTION: PROCESS"):
            current_section = "PROCESS"
        elif line.startswith("PATH_START"):
            current_path = []
        elif line.startswith("END"):
            if current_path:
                paths.append(current_path)
            current_path = []

        # --- 解析结果线 (交线 OR 最近点) ---
        elif line.startswith("RESULT_INTERSECTION") or line.startswith("RESULT_NEAREST"):
            parts = line.split()
            # 格式: RESULT_XXX Ax Ay Az Bx By Bz Dist
            res_type = 'intersection' if "INTERSECTION" in line else 'nearest'
            res = {
                'type': 'result',
                'res_type': res_type,
                'pA': [float(parts[1]), float(parts[2]), float(parts[3])],
                'pB': [float(parts[4]), float(parts[5]), float(parts[6])],
                'dist': float(parts[7])
            }
            current_path.append(res)

        # --- 解析搜索步骤 ---
        elif line.startswith("STEP"):
            parts = line.split()
            if len(parts) < 11: continue
            step = {
                'type': 'step',
                'A': {'type': int(parts[2]), 'u': int(parts[3]), 'v': int(parts[4])},
                'B': {'type': int(parts[6]), 'u': int(parts[7]), 'v': int(parts[8])},
                'dist': float(parts[10])
            }
            current_path.append(step)

        # --- 解析几何点数据 ---
        else:
            if current_section in ["A", "B"]:
                try:
                    coords = list(map(float, line.split()))
                    temp_points.append(coords)
                    if dims[0] > 0 and len(temp_points) == dims[0] * dims[1]:
                        grid = np.array(temp_points).reshape(dims[1], dims[0], 3)
                        if current_section == "A":
                            data_A = grid
                        else:
                            data_B = grid
                except ValueError:
                    continue

    return data_A, data_B, paths


# ==========================================
# 2. 几何辅助函数 (修正后的 Type 2 逻辑)
# ==========================================
def get_tri_coords(grid, t_type, u, v):
    rows, cols, _ = grid.shape
    # 边界保护
    u = max(0, min(u, cols - 1))
    v = max(0, min(v, rows - 1))
    u_next = min(u + 1, cols - 1)
    v_next = min(v + 1, rows - 1)

    if t_type == 1:
        # 左下: (u,v) -> (u+1,v) -> (u,v+1)
        p1 = grid[v, u]
        p2 = grid[v, u_next]
        p3 = grid[v_next, u]
    else:
        # 右上: (u+1,v) -> (u+1,v+1) -> (u,v+1)
        p1 = grid[v, u_next]
        p2 = grid[v_next, u_next]
        p3 = grid[v_next, u]

    # 返回闭合坐标
    return [p1[0], p2[0], p3[0], p1[0]], \
        [p1[1], p2[1], p3[1], p1[1]], \
        [p1[2], p2[2], p3[2], p1[2]]


# 辅助：创建线框网格
def create_mesh_wireframe(grid, color, name, opacity=0.3):
    x, y, z = [], [], []
    rows, cols, _ = grid.shape
    for r in range(rows):
        for c in range(cols):
            x.append(grid[r, c, 0]);
            y.append(grid[r, c, 1]);
            z.append(grid[r, c, 2])
        x.append(None);
        y.append(None);
        z.append(None)
    for c in range(cols):
        for r in range(rows):
            x.append(grid[r, c, 0]);
            y.append(grid[r, c, 1]);
            z.append(grid[r, c, 2])
        x.append(None);
        y.append(None);
        z.append(None)
    return go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color=color, width=1), opacity=opacity, name=name)


# 辅助：创建实体曲面
def create_mesh_surface(grid, color, name, opacity=0.8):
    return go.Surface(
        x=grid[:, :, 0], y=grid[:, :, 1], z=grid[:, :, 2],
        opacity=opacity, name=name, showscale=False,
        colorscale=[[0, color], [1, color]]
    )


# ==========================================
# 3. 主程序
# ==========================================
def main():
    # 1. 设置路径
    filename = "Main/Main/debug_output.txt"

    print(f"Start parsing: {filename}")
    data_A, data_B, paths = parse_file(filename)

    if data_A is None or not paths:
        print("Data load failed or no paths found.")
        return

    # ---------------------------------------------------------
    # 准备工作：提取所有结果线段用于 Summary
    # ---------------------------------------------------------
    inter_x, inter_y, inter_z = [], [], []  # 交线 (红)
    near_x, near_y, near_z = [], [], []  # 最近线 (金)

    for path in paths:
        # 查找该路径中是否包含结果行
        res = next((item for item in path if item['type'] == 'result'), None)
        if res:
            if res['res_type'] == 'intersection':
                inter_x.extend([res['pA'][0], res['pB'][0], None])
                inter_y.extend([res['pA'][1], res['pB'][1], None])
                inter_z.extend([res['pA'][2], res['pB'][2], None])
            else:
                near_x.extend([res['pA'][0], res['pB'][0], None])
                near_y.extend([res['pA'][1], res['pB'][1], None])
                near_z.extend([res['pA'][2], res['pB'][2], None])

    trace_inter = go.Scatter3d(
        x=inter_x, y=inter_y, z=inter_z, mode='lines',
        line=dict(color='red', width=8), name='Intersections'
    )
    trace_near = go.Scatter3d(
        x=near_x, y=near_y, z=near_z, mode='lines',
        line=dict(color='gold', width=4), name='Nearest Links'
    )

    # ---------------------------------------------------------
    # 输出 1: 曲面模式总结 (summary_surface.html)
    # ---------------------------------------------------------
    surf_A = create_mesh_surface(data_A, 'blue', 'Surface A', 0.6)
    surf_B = create_mesh_surface(data_B, 'red', 'Surface B', 0.6)

    fig_surf = go.Figure(data=[surf_A, surf_B, trace_inter, trace_near])
    fig_surf.update_layout(
        title="Summary: Solid Surfaces + Results (Red=Intersect, Gold=Nearest)",
        scene=dict(aspectmode='data')
    )
    fig_surf.write_html("summary_surface.html")
    print("Generated: summary_surface.html")

    # ---------------------------------------------------------
    # 输出 2: 线框模式总结 (summary_wireframe.html)
    # ---------------------------------------------------------
    wire_A = create_mesh_wireframe(data_A, 'blue', 'Grid A', 0.5)
    wire_B = create_mesh_wireframe(data_B, 'red', 'Grid B', 0.5)

    fig_wire = go.Figure(data=[wire_A, wire_B, trace_inter, trace_near])
    fig_wire.update_layout(
        title="Summary: Wireframes + Results (Red=Intersect, Gold=Nearest)",
        scene=dict(aspectmode='data')
    )
    fig_wire.write_html("summary_wireframe.html")
    print("Generated: summary_wireframe.html")

    # ---------------------------------------------------------
    # 输出 3: 详细动画 (greedy_viz_animation.html)
    # ---------------------------------------------------------
    print("Generating Animation Frames (this may take a moment)...")

    # 基础背景
    anim_bg_A = create_mesh_wireframe(data_A, 'blue', 'Mesh A', 0.2)
    anim_bg_B = create_mesh_wireframe(data_B, 'red', 'Mesh B', 0.2)

    # 占位符 Trace
    # 历史路径 (灰色)
    trace_hist_A = go.Scatter3d(x=[], y=[], z=[], mode='lines', line=dict(color='#888888', width=2), name='Hist A')
    trace_hist_B = go.Scatter3d(x=[], y=[], z=[], mode='lines', line=dict(color='#888888', width=2), name='Hist B')
    # 当前三角片 (高亮)
    trace_curr_A = go.Scatter3d(x=[], y=[], z=[], mode='lines', line=dict(color='cyan', width=5), name='Curr A')
    trace_curr_B = go.Scatter3d(x=[], y=[], z=[], mode='lines', line=dict(color='magenta', width=5), name='Curr B')
    # 搜索连线 (绿色虚线)
    trace_search_link = go.Scatter3d(x=[], y=[], z=[], mode='lines', line=dict(color='green', width=4, dash='dash'),
                                     name='Search Dist')
    # 最终结果线 (每条路径单独显示，金色或红色)
    trace_path_result = go.Scatter3d(x=[], y=[], z=[], mode='lines', line=dict(width=10), name='Path Result')

    frames = []
    sliders_steps = []
    path_buttons = []

    for p_idx, path in enumerate(paths):
        hist_Ax, hist_Ay, hist_Az = [], [], []
        hist_Bx, hist_By, hist_Bz = [], [], []

        path_frame_names = []

        # 获取该路径的最终结果（如果有）
        path_res = next((item for item in path if item['type'] == 'result'), None)

        # 遍历步骤
        steps_only = [item for item in path if item['type'] == 'step']

        for s_idx, step in enumerate(steps_only):
            # 1. 计算当前三角形坐标
            ax, ay, az = get_tri_coords(data_A, step['A']['type'], step['A']['u'], step['A']['v'])
            bx, by, bz = get_tri_coords(data_B, step['B']['type'], step['B']['u'], step['B']['v'])

            # 2. 计算中心点 (用于绿色连线)
            acx, acy, acz = np.mean(ax[:3]), np.mean(ay[:3]), np.mean(az[:3])
            bcx, bcy, bcz = np.mean(bx[:3]), np.mean(by[:3]), np.mean(bz[:3])

            # 3. 更新历史
            hist_Ax.extend(ax + [None])
            hist_Ay.extend(ay + [None])
            hist_Az.extend(az + [None])
            hist_Bx.extend(bx + [None])
            hist_By.extend(by + [None])
            hist_Bz.extend(bz + [None])

            # 4. 准备结果线 (仅在最后一步显示)
            res_x, res_y, res_z = [], [], []
            res_color = 'gold'
            is_last_step = (s_idx == len(steps_only) - 1)

            if is_last_step and path_res:
                res_x = [path_res['pA'][0], path_res['pB'][0]]
                res_y = [path_res['pA'][1], path_res['pB'][1]]
                res_z = [path_res['pA'][2], path_res['pB'][2]]
                res_color = 'red' if path_res['res_type'] == 'intersection' else 'gold'

            frame_name = f"P{p_idx}_S{s_idx}"
            path_frame_names.append(frame_name)

            frames.append(go.Frame(
                name=frame_name,
                data=[
                    anim_bg_A, anim_bg_B,
                    go.Scatter3d(x=hist_Ax, y=hist_Ay, z=hist_Az),  # Hist A
                    go.Scatter3d(x=hist_Bx, y=hist_By, z=hist_Bz),  # Hist B
                    go.Scatter3d(x=ax, y=ay, z=az),  # Curr A
                    go.Scatter3d(x=bx, y=by, z=bz),  # Curr B
                    go.Scatter3d(x=[acx, bcx], y=[acy, bcy], z=[acz, bcz]),  # Search Link
                    go.Scatter3d(x=res_x, y=res_y, z=res_z, line=dict(color=res_color))  # Result Line
                ],
                layout=go.Layout(
                    title_text=f"Path {p_idx} | Step {s_idx} | Dist: {step['dist']:.4f}"
                )
            ))

            sliders_steps.append({
                "args": [[frame_name], {"frame": {"duration": 100, "redraw": True}, "mode": "immediate"}],
                "label": f"{p_idx}-{s_idx}", "method": "animate"
            })

        # 下拉菜单按钮
        path_buttons.append({
            "args": [path_frame_names, {"frame": {"duration": 300, "redraw": True}, "mode": "immediate"}],
            "label": f"Path {p_idx}", "method": "animate"
        })

    # 动画 Layout
    initial_data = [anim_bg_A, anim_bg_B, trace_hist_A, trace_hist_B, trace_curr_A, trace_curr_B, trace_search_link,
                    trace_path_result]

    layout_anim = go.Layout(
        title="Greedy Search Animation",
        scene=dict(xaxis=dict(showgrid=False), yaxis=dict(showgrid=False), zaxis=dict(showgrid=False),
                   aspectmode='data'),
        updatemenus=[
            {"buttons": path_buttons, "direction": "down", "pad": {"r": 10, "t": 10}, "x": 0.1, "y": 1.1,
             "xanchor": "left", "yanchor": "top"},
            {"buttons": [
                {"args": [None, {"frame": {"duration": 100, "redraw": True}, "fromcurrent": True}], "label": "Play",
                 "method": "animate"},
                {"args": [[None], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}], "label": "Pause",
                 "method": "animate"}
            ],
                "direction": "left", "pad": {"r": 10, "t": 87}, "x": 0.1, "y": 0, "xanchor": "right", "yanchor": "top"}
        ],
        sliders=[{"active": 0, "steps": sliders_steps, "x": 0.1, "y": 0}]
    )

    fig_anim = go.Figure(data=initial_data, layout=layout_anim, frames=frames)
    fig_anim.write_html("greedy_viz_animation.html")
    print("Generated: greedy_viz_animation.html")
    print("All tasks completed.")


if __name__ == "__main__":
    main()