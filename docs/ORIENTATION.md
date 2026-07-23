# Orientation Handling

This document explains the orientation parameters used by `image.sh`,
`filter.sh`, and clone analysis, and how to choose the correct setting with the
schematic images.

Pass both required values to the image step:

- `--orientation <mode>`: one of `normal`, `horizontal`, `vertical`, or `rotate`
- `--swap-xy <bool>`: `True` swaps the x/y axes after applying orientation

These parameters are used when the image coordinate system, spatial barcode coordinate system, and final visualization direction do not match.

`dbit.sh image` appends both values to the per-dataset config. The later filter
and clone steps consume those stored values and do not accept separate
orientation arguments.

## 1. Where Orientation Is Applied

### `image.sh`

`image.sh` splits `align.png` according to the spatial grid and names each tile as `x_y.tif`. The orientation parameters determine where spots such as `0_0`, `1_0`, and `0_1` are mapped on the original image.

For example:

```bash
bash dbit.sh image --input /path/to/align.png --config /path/to/dbit.config.sh \
    --orientation horizontal --swap-xy False
```

### `filter.sh`

`filter.sh` merges transcriptome/amplicon filtered plots onto `gray.png`. The orientation parameters determine how the filtered plots are flipped or rotated before merging, so that they align with `gray.png`.

### Clone analysis

Clone analysis reads the same stored `orientation` and `swap_xy` values. It
applies them to the mRNA Leiden-cluster spatial background and clone spot
coordinates so the top-LR grid has the same direction as the image and plot
steps.

Clone analysis also provides a separate presentation option:

```bash
dbit clone --config /path/to/dbit.config.sh --rotate 90
```

`--rotate` accepts `0`, `90`, `180`, or `270` degrees clockwise. This option
rotates only the final spatial UMAP/Leiden grid drawn by `top_lr_plot.py`. It
does not rotate the photographed tissue image, change `gray.png`, or rotate the
plot legends. The value is stored as `clone_rotate` in the dataset config.

Use `orientation` and `swap_xy` to align data coordinates with the photographed
image. Use clone `--rotate` only when the already aligned grid needs a different
direction for presentation.

## 2. Orientation Options

### `normal`

No orientation adjustment is applied. This assumes that the spatial coordinates and image direction already match.

Use this when:

- spot `0_0` is already in the expected position;
- the x and y coordinate directions match the displayed image;
- filtered plots align with `gray.png` without any transformation.

Command:

```bash
orientation=normal
swap_xy=False
```

### `horizontal`

Flip the coordinate system or overlay image horizontally.

![Horizontal orientation](image/horizontal.png)

Use this when:

- the left-right direction of the image is reversed relative to the spatial coordinates;
- spot `0_0` is on the opposite side in the horizontal direction;
- the overlay and `gray.png` show a left-right mirror mismatch.

Command:

```bash
orientation=horizontal
swap_xy=False
```

### `vertical`

Flip the coordinate system or overlay image vertically.

![Vertical orientation](image/vertical.png)

Use this when:

- the top-bottom direction of the image is reversed relative to the spatial coordinates;
- spot `0_0` is on the opposite side in the vertical direction;
- the overlay and `gray.png` show a top-bottom mirror mismatch.

Command:

```bash
orientation=vertical
swap_xy=False
```

### `rotate`

Rotate by 180 degrees. This is equivalent to applying both `horizontal` and `vertical`.

![Rotate orientation](image/rotate.png)

Use this when:

- both the left-right and top-bottom directions are reversed;
- spot `0_0` is located at the opposite corner;
- the overlay and `gray.png` differ by 180 degrees.

Command:

```bash
orientation=rotate
swap_xy=False
```

## 3. 90-Degree Rotation

90-degree rotation is handled by combining `orientation` with `swap_xy=True`.

### 90 Degrees Counterclockwise

Use `horizontal` with `swap_xy=True`.

![Counterclockwise orientation](image/counterclockwise.png)

Command:

```bash
orientation=horizontal
swap_xy=True
```

Use this when:

- the image needs to be rotated 90 degrees counterclockwise relative to the spatial coordinates;
- the x/y axes are swapped and the origin position needs to be corrected by a horizontal flip.

### 90 Degrees Clockwise

Use `vertical` with `swap_xy=True`.

![Clockwise orientation](image/clockwise.png)

Command:

```bash
orientation=vertical
swap_xy=True
```

Use this when:

- the image needs to be rotated 90 degrees clockwise relative to the spatial coordinates;
- the x/y axes are swapped and the origin position needs to be corrected by a vertical flip.

## 4. Quick Reference

| Required transformation | Parameters |
| --- | --- |
| No adjustment | `orientation=normal`, `swap_xy=False` |
| Flip left-right | `orientation=horizontal`, `swap_xy=False` |
| Flip top-bottom | `orientation=vertical`, `swap_xy=False` |
| Rotate 180 degrees | `orientation=rotate`, `swap_xy=False` |
| Rotate 90 degrees counterclockwise | `orientation=horizontal`, `swap_xy=True` |
| Rotate 90 degrees clockwise | `orientation=vertical`, `swap_xy=True` |

Clone presentation rotation is independent of this table:

| Presentation result | Clone option |
| --- | --- |
| Keep the aligned grid direction | `--rotate 0` |
| Rotate the grid 90 degrees clockwise | `--rotate 90` |
| Rotate the grid 180 degrees | `--rotate 180` |
| Rotate the grid 270 degrees clockwise | `--rotate 270` |

## 5. Recommended Check

1. Run `image.sh` first and inspect the generated `result.png`.
2. Check whether spot `0_0` in `result.png` is in the expected position.
3. If the labels are incorrect, rerun the image step with adjusted `--orientation` and `--swap-xy` values.
4. Run `filter.sh` and check whether the generated `merged_*_filtered.png` files align with `gray.png`.
5. Run clone analysis and confirm that the clone spots and Leiden background
   have the same orientation. If only the presentation direction needs to
   change, rerun clone with `--rotate`; do not change the image orientation.

Use the same orientation parameters for `image.sh`, `filter.sh`, and clone
analysis. This keeps image splitting, cell filtering, merged plots, and clone
coordinates consistent.
