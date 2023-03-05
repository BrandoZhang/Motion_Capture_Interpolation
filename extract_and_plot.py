import os
import pyecharts.options as opts
from pyecharts.charts import Line


def extract_data(amc_filename):
    lfemur_x_series = []
    root_z_series = []

    with open(amc_filename, 'r') as f:
        for line in f:
            if line.startswith("lfemur"):
                lfemur_x_series.append(float(line.split()[1]))  # Add X rotation
            elif line.startswith("root"):
                root_z_series.append(float(line.split()[6]))  # Add Z rotation

    return lfemur_x_series, root_z_series


def render_HTML_graph(title,
                      subtitle,
                      save_path,
                      method1,
                      method1_data,
                      method2,
                      method2_data,
                      method3,
                      method3_data,
                      frame_range):

    (
        Line()
        .add_xaxis(xaxis_data=[str(i) for i in frame_range])
        .add_yaxis(
            series_name=method1,
            y_axis=method1_data[frame_range.start: frame_range.stop],
            label_opts=opts.LabelOpts(is_show=False),
        )
        .add_yaxis(
            series_name=method2,
            y_axis=method2_data[frame_range.start: frame_range.stop],
            label_opts=opts.LabelOpts(is_show=False),
        )
        .add_yaxis(
            series_name=method3,
            y_axis=method3_data[frame_range.start: frame_range.stop],
            label_opts=opts.LabelOpts(is_show=False),
        )
        .set_global_opts(
            xaxis_opts=opts.AxisOpts(name="Frame"),
            # yaxis_opts=opts.AxisOpts(name="Euler angle"),
            title_opts=opts.TitleOpts(title=title, subtitle=subtitle),
            tooltip_opts=opts.TooltipOpts(trigger="axis"),
            toolbox_opts=opts.ToolboxOpts(is_show=True),
        )
        .render(save_path)
    )


input_amc = os.path.join('.', 'docs/graph', '131_04-dance.amc')
linear_euler_amc = os.path.join('.', 'docs/graph', 'outputLE20DanceMotion.amc')
bezier_euler_amc = os.path.join('.', 'docs/graph', 'outputBE20DanceMotion.amc')
slerp_amc = os.path.join('.', 'docs/graph', 'outputLQ20DanceMotion.amc')
bezier_slerp_amc = os.path.join('.', 'docs/graph', 'outputBQ20DanceMotion.amc')

lfemur_frame_range = range(600, 800)
root_frame_range = range(200, 500)

input_lfemur_x_series, input_root_z_series = extract_data(input_amc)
linear_euler_lfemur_x_series, linear_euler_root_z_series = extract_data(linear_euler_amc)
bezier_euler_lfemur_x_series, bezier_euler_root_z_series = extract_data(bezier_euler_amc)
slerp_lfemur_x_series, slerp_root_z_series = extract_data(slerp_amc)
bezier_slerp_lfemur_x_series, bezier_slerp_root_z_series = extract_data(bezier_slerp_amc)

# Graph 1: Linear Euler vs Bezier Euler vs Input
# `131_04-dance.amc` `lfemur` joint around X axis, frames 600-800, N=20
render_HTML_graph(title="Graph 1",
                  subtitle="Linear Euler vs Bezier Euler vs Input",
                  save_path=os.path.join('.', 'docs/graph', 'graph1.html'),
                  method1="Linear Euler",
                  method1_data=linear_euler_lfemur_x_series,
                  method2="Bezier Euler",
                  method2_data=bezier_euler_lfemur_x_series,
                  method3="Input",
                  method3_data=input_lfemur_x_series,
                  frame_range=lfemur_frame_range)


# Graph 2: SLERP vs Bezier SLERP vs Input
# `131_04-dance.amc` `lfemur` joint around X axis, frames 600-800, N=20
render_HTML_graph(title="Graph 2",
                  subtitle="SLERP vs Bezier SLERP vs Input",
                  save_path=os.path.join('.', 'docs/graph', 'graph2.html'),
                  method1="SLERP",
                  method1_data=slerp_lfemur_x_series,
                  method2="Bezier SLERP",
                  method2_data=bezier_slerp_lfemur_x_series,
                  method3="Input",
                  method3_data=input_lfemur_x_series,
                  frame_range=lfemur_frame_range)

# Graph 3: Linear Euler vs SLERP vs Input
# `131_04-dance.amc` `root` joint around Z axis, frames 200-500, N=20
render_HTML_graph(title="Graph 3",
                  subtitle="Linear Euler vs SLERP vs Input",
                  save_path=os.path.join('.', 'docs/graph', 'graph3.html'),
                  method1="Linear Euler",
                  method1_data=linear_euler_root_z_series,
                  method2="SLERP",
                  method2_data=slerp_root_z_series,
                  method3="Input",
                  method3_data=input_root_z_series,
                  frame_range=root_frame_range)

# Graph 4: Bezier Euler vs Bezier SLERP vs Input
# `131_04-dance.amc` `root` joint around Z axis, frames 200-500, N=20
render_HTML_graph(title="Graph 4",
                  subtitle="Bezier Euler vs Bezier SLERP vs Input",
                  save_path=os.path.join('.', 'docs/graph', 'graph4.html'),
                  method1="Bezier Euler",
                  method1_data=bezier_euler_root_z_series,
                  method2="Bezier SLERP",
                  method2_data=bezier_slerp_root_z_series,
                  method3="Input",
                  method3_data=input_root_z_series,
                  frame_range=root_frame_range)
