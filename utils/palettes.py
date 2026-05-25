"""Generate colorblind-safe color palettes for cell types and brain regions."""

import glasbey as gb
import seaborn as sns

# Muted palette for cell types (suitable for printed figures)
muted = gb.create_block_palette(
    [4, 5, 2, 2, 3, 1, 1, 1, 1, 1, 1],
    lightness_bounds=(10, 70),
    chroma_bounds=(0, 70),
    hue_bounds=(0, 360),
    colorblind_safe=True,
    optimize_palette=True,
)
sns.palplot(muted)

# Vibrant palette for cell types (suitable for digital figures)
vibrant = gb.create_block_palette(
    [4, 5, 2, 2, 3, 1, 1, 1, 1, 1, 1],
    lightness_bounds=(0, 100),
    chroma_bounds=(0, 150),
    hue_bounds=(0, 360),
    colorblind_safe=True,
    optimize_palette=True,
)
sns.palplot(vibrant)

# Brain region palette (5 regions)
sns.palplot(
    gb.create_palette(
        5,
        lightness_bounds=(30, 80),
        chroma_bounds=(0, 55),
        hue_bounds=(0, 360),
        colorblind_safe=True,
        optimize_palette=True,
    )
)
gb.create_palette(
    5,
    lightness_bounds=(30, 80),
    chroma_bounds=(0, 55),
    hue_bounds=(0, 360),
    colorblind_safe=True,
    optimize_palette=True,
)
