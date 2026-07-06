import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from PIL import Image

def legend(data, output_path, name):
    values = np.asarray(data, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if values.size == 0:
        return
    p95 = np.percentile(values, 95)

    red_rgba = np.zeros((256, 4))
    red_rgba[:, 0] = 1.0  # 红色通道
    red_rgba[:, 1] = 0.0  # 绿色
    red_rgba[:, 2] = 0.0  # 蓝色
    red_rgba[:, 3] = np.linspace(0, 1, 256)  # 透明度从 0 渐变到 1
    custom_red_alpha_cmap = mcolors.ListedColormap(red_rgba)

    norm = mcolors.Normalize(vmin=0, vmax=p95)
    sm = plt.cm.ScalarMappable(cmap=custom_red_alpha_cmap, norm=norm)
    sm.set_array([])
    fig = plt.figure(figsize=(5, 10))
    ax = fig.add_axes([0.10, 0.14, 0.14, 0.70])
    cbar = fig.colorbar(sm, cax=ax, orientation='vertical')
    ticks = np.linspace(0, p95, 5)
    cbar.set_ticks(ticks)
    if p95 >= 10:
        ticklabels = [f"{t:,.0f}" for t in ticks]
        p95_label = f"{p95:,.0f}"
    else:
        ticklabels = [f"{t:.2f}" for t in ticks]
        p95_label = f"{p95:.2f}"
    cbar.set_ticklabels(ticklabels, fontsize=22)
    cbar.ax.tick_params(labelsize=22, pad=10)
    display_name = name.replace('_', ' ')
    cbar.set_label(display_name, fontsize=22, labelpad=28)
    fig.suptitle(
        f"Color capped at P95\nP95 = {p95_label}",
        fontsize=22,
        y=0.96,
    )
    fig.text(
        0.5,
        0.04,
        "Values ≥ P95 use maximum opacity",
        ha='center',
        va='center',
        fontsize=16,
    )
    ax.tick_params(labelsize=22)
    fig.savefig(f'{output_path}/{name}', dpi=300)
    plt.close()

def plot_frame(data, output_path, config):

    positive_cell_spots = data['count'] > 0 if 'count' in data.columns else None
    has_umi_per_cell = (
        'umi_count' in data.columns
        and positive_cell_spots is not None
        and bool((positive_cell_spots & (data['umi_count'] > 0)).any())
    )
    has_gene_per_cell = (
        'gene_count' in data.columns
        and positive_cell_spots is not None
        and bool((positive_cell_spots & (data['gene_count'] > 0)).any())
    )

    if 'umi_count' in data.columns:
        legend(data['umi_count'], output_path, 'UMI_distribution')
        if has_umi_per_cell:
            legend(data.loc[positive_cell_spots, 'umi_count'] / data.loc[positive_cell_spots, 'count'], output_path, 'UMI_per_cell_distribution')
    if 'gene_count' in data.columns:
        legend(data['gene_count'], output_path, 'Gene_distribution')
        if has_gene_per_cell:
            legend(data.loc[positive_cell_spots, 'gene_count'] / data.loc[positive_cell_spots, 'count'], output_path, 'Gene_per_cell_distribution')

    if 'leiden' in data.columns:
        frame = np.zeros((int(config.x_spots_number * config.length_spot + (config.x_spots_number-1) * config.interval),
                          int(config.y_spots_number * config.length_spot + (config.y_spots_number-1) * config.interval), 4),
                          dtype=np.uint8)
        if 'color' in data.columns:
            cluster_colors = {
                int(row['leiden']): np.array(mcolors.to_rgba(row['color']))
                for _, row in data[['leiden', 'color']].drop_duplicates().iterrows()
            }
        else:
            clusters = sorted(data['leiden'].astype(int).unique())
            cmap = plt.get_cmap('tab20', len(clusters))
            cluster_colors = {
                cluster_id: np.array(cmap(i))
                for i, cluster_id in enumerate(clusters)
            }

    if 'umi_count' in data.columns:
        frame_umi = np.zeros((int(config.x_spots_number * config.length_spot + (config.x_spots_number-1) * config.interval),
                              int(config.y_spots_number * config.length_spot + (config.y_spots_number-1) * config.interval), 4),
                              dtype=np.float32)
        if has_umi_per_cell:
            frame_umi_per_cell = np.zeros((int(config.x_spots_number * config.length_spot + (config.x_spots_number-1) * config.interval),
                                           int(config.y_spots_number * config.length_spot + (config.y_spots_number-1) * config.interval), 4),
                                           dtype=np.float32)

    if 'gene_count' in data.columns:
        frame_gene = np.zeros((int(config.x_spots_number * config.length_spot + (config.x_spots_number-1) * config.interval),
                               int(config.y_spots_number * config.length_spot + (config.y_spots_number-1) * config.interval), 4),
                               dtype=np.float32)
        if has_gene_per_cell:
            frame_gene_per_cell = np.zeros((int(config.x_spots_number * config.length_spot + (config.x_spots_number-1) * config.interval),
                                            int(config.y_spots_number * config.length_spot + (config.y_spots_number-1) * config.interval), 4),
                                            dtype=np.float32)

    for _, row in data.iterrows():
        x_idx = int(row['x'])
        y_idx = int(row['y'])

        x_start = x_idx * (config.length_spot + config.interval)
        y_start = y_idx * (config.length_spot + config.interval)
        x_end = x_start + config.length_spot
        y_end = y_start + config.length_spot

        if 'umi_count' in data.columns:
            frame_umi[y_start: y_end, x_start: x_end, 0] = 255
            frame_umi[y_start: y_end, x_start: x_end, 3] = row['umi_count']
            if has_umi_per_cell and row['count'] > 0:
                frame_umi_per_cell[y_start: y_end, x_start: x_end, 0] = 255
                frame_umi_per_cell[y_start: y_end, x_start: x_end, 3] = row['umi_count'] / row['count']
        if 'gene_count' in data.columns:
            frame_gene[y_start: y_end, x_start: x_end, 0] = 255
            frame_gene[y_start: y_end, x_start: x_end, 3] = int(row['gene_count'])
            if has_gene_per_cell and row['count'] > 0:
                frame_gene_per_cell[y_start: y_end, x_start: x_end, 0] = 255
                frame_gene_per_cell[y_start: y_end, x_start: x_end, 3] = int(row['gene_count']) / row['count']
        if 'leiden' in data.columns:
            cluster_id = int(row['leiden'])
            frame[y_start: y_end, x_start: x_end, :] = (cluster_colors[cluster_id] * 255).astype(np.uint8)

    if 'leiden' in data.columns:
        img_umap = Image.fromarray(frame, mode = 'RGBA')
        img_umap = img_umap.resize((int((config.x_spots_number * config.length_spot + (config.x_spots_number - 1) * config.interval)/config.pixel_length),
                                    int((config.y_spots_number * config.length_spot + (config.y_spots_number - 1) * config.interval)/config.pixel_length)),
                                    resample = Image.NEAREST)
        img_umap.save(f'{output_path}/umap_filtered.png')
    if 'umi_count' in data.columns:
        nozero = frame_umi[:, :, 3][frame_umi[:, :, 3] > 0]
        medain = np.percentile(nozero, 95)
        frame_umi[:, :, 3] = frame_umi[:, :, 3] / medain
        frame_umi[:, :, 3][frame_umi[:, :, 3] > 1] = 1
        frame_umi[:, :, 3] = (frame_umi[:, :, 3] * 255).astype(np.uint8)
        frame_umi = frame_umi.astype(np.uint8)

        if has_umi_per_cell:
            nozero = frame_umi_per_cell[:, :, 3][frame_umi_per_cell[:, :, 3] > 0]
            medain = np.percentile(nozero, 95)
            frame_umi_per_cell[:, :, 3] = frame_umi_per_cell[:, :, 3] / medain
            frame_umi_per_cell[:, :, 3][frame_umi_per_cell[:, :, 3] > 1] = 1
            frame_umi_per_cell[:, :, 3] = (frame_umi_per_cell[:, :, 3] * 255).astype(np.uint8)
            frame_umi_per_cell = frame_umi_per_cell.astype(np.uint8)

        img_umi = Image.fromarray(frame_umi, mode = 'RGBA')
        img_umi = img_umi.resize((int((config.x_spots_number * config.length_spot + (config.x_spots_number - 1) * config.interval)/config.pixel_length),
                                  int((config.y_spots_number * config.length_spot + (config.y_spots_number - 1) * config.interval)/config.pixel_length)),
                                  resample = Image.NEAREST)
        img_umi.save(f'{output_path}/umi_filtered.png')

        if has_umi_per_cell:
            img_umi_per_cell = Image.fromarray(frame_umi_per_cell, mode = 'RGBA')
            img_umi_per_cell = img_umi_per_cell.resize((int((config.x_spots_number * config.length_spot + (config.x_spots_number - 1) * config.interval)/config.pixel_length),
                                                        int((config.y_spots_number * config.length_spot + (config.y_spots_number - 1) * config.interval)/config.pixel_length)),
                                                        resample = Image.NEAREST)
            img_umi_per_cell.save(f'{output_path}/umi_per_cell_filtered.png')

    if 'gene_count' in data.columns:
        nozero = frame_gene[:, :, 3][frame_gene[:, :, 3] > 0]
        medain = np.percentile(nozero, 95)
        frame_gene[:, :, 3] = frame_gene[:, :, 3] / medain
        frame_gene[:, :, 3][frame_gene[:, :, 3] > 1] = 1
        frame_gene[:, :, 3] = (frame_gene[:, :, 3] * 255).astype(np.uint8)
        frame_gene = frame_gene.astype(np.uint8)

        if has_gene_per_cell:
            nozero = frame_gene_per_cell[:, :, 3][frame_gene_per_cell[:, :, 3] > 0]
            medain = np.percentile(nozero, 95)
            frame_gene_per_cell[:, :, 3] = frame_gene_per_cell[:, :, 3] / medain
            frame_gene_per_cell[:, :, 3][frame_gene_per_cell[:, :, 3] > 1] = 1
            frame_gene_per_cell[:, :, 3] = (frame_gene_per_cell[:, :, 3] * 255).astype(np.uint8)
            frame_gene_per_cell = frame_gene_per_cell.astype(np.uint8)

        
        img_gene = Image.fromarray(frame_gene, mode = 'RGBA')
        img_gene = img_gene.resize((int((config.x_spots_number * config.length_spot + (config.x_spots_number - 1) * config.interval)/config.pixel_length),
                                    int((config.y_spots_number * config.length_spot + (config.y_spots_number - 1) * config.interval)/config.pixel_length)),
                                    resample = Image.NEAREST)
        img_gene.save(f'{output_path}/gene_filtered.png')

        if has_gene_per_cell:
            img_gene_per_cell = Image.fromarray(frame_gene_per_cell, mode = 'RGBA')
            img_gene_per_cell = img_gene_per_cell.resize((int((config.x_spots_number * config.length_spot + (config.x_spots_number - 1) * config.interval)/config.pixel_length),
                                                          int((config.y_spots_number * config.length_spot + (config.y_spots_number - 1) * config.interval)/config.pixel_length)),
                                                          resample = Image.NEAREST)
            img_gene_per_cell.save(f'{output_path}/gene_per_cell_filtered.png')
