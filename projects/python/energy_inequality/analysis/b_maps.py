"""Map helpers shared by the b_ scripts (choropleths of clubs and local
coefficients, with a London inset so small urban units stay visible).

Figures are 5.4 in wide so that 12-13 pt type renders at roughly 12-13 px
when the image is scaled to a 375 px phone column."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Rectangle

import b_common as bc

# Bounding box of Greater London (BNG metres) for the inset.
LONDON = (503000, 155000, 562000, 201000)
HATCH_VALUE = 0  # category drawn hatched (divergent units)


def _draw(ax, gdf, column, colours=None, cmap=None, norm=None, ns_mask=None, edge=0.0):
    if colours is not None:
        fill = gdf[column].map(colours).fillna("#ffffff")
        hatched = gdf[column] == HATCH_VALUE
        gdf[~hatched].plot(ax=ax, color=fill[~hatched], linewidth=edge, edgecolor="white")
        if hatched.any():
            gdf[hatched].plot(ax=ax, facecolor="white", edgecolor=bc.INK, hatch="////", linewidth=0.3)
    else:
        mask = None if ns_mask is None else ns_mask.reindex(gdf.index).fillna(False).values.astype(bool)
        if mask is not None and mask.any():
            gdf[mask].plot(ax=ax, color=bc.NS_GREY, linewidth=0)
        (gdf if mask is None else gdf[~mask]).plot(ax=ax, column=column, cmap=cmap, norm=norm,
                                                   linewidth=edge, edgecolor="white")
    ax.set_axis_off()


def choropleth(gdf, column, title, name, colours=None, labels=None, cmap=None, norm=None,
               ns_mask=None, cbar_label=None, boundaries=None, inset=True, note=None, cbar_ticks=None):
    """Map of England and Wales with an optional London inset.

    colours: dict value -> colour for categorical/ordinal maps (legend uses
    `labels`; value 0 is drawn hatched); otherwise a continuous `cmap`/`norm`
    with a colourbar. ns_mask: boolean Series (True = not significant)
    drawn in light grey."""
    n_leg = len(colours) + (ns_mask is not None) if colours is not None else 0
    leg_h = 0.0 if colours is None else 0.045 * int(np.ceil(n_leg / 3))
    fig_h = 7.0 + leg_h * 7.0
    fig = plt.figure(figsize=(5.4, fig_h))
    bottom = (0.20 if colours is None else leg_h + 0.06) * 7.0 / fig_h
    ax = fig.add_axes([0.0, bottom, 1.0, 1 - bottom - 0.10 * 7.0 / fig_h])
    edge = 0.15 if len(gdf) < 1000 else 0.0
    _draw(ax, gdf, column, colours, cmap, norm, ns_mask, edge)
    if boundaries is not None:
        boundaries.boundary.plot(ax=ax, color="white", linewidth=0.3)
    fig.text(0.01, 1 - 0.015 * 7.0 / fig_h, title, fontsize=14, va="top", ha="left")
    if inset:
        x0, y0, x1, y1 = LONDON
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=bc.INK, lw=0.8))
        pos = ax.get_position()
        axi = fig.add_axes([0.63, pos.y0 + pos.height * 0.63, 0.36, pos.height * 0.27])
        sub = gdf.cx[x0:x1, y0:y1]
        _draw(axi, sub, column, colours, cmap, norm,
              None if ns_mask is None else ns_mask.reindex(sub.index), 0.1 if len(gdf) < 9000 else 0.0)
        if boundaries is not None:
            boundaries.cx[x0:x1, y0:y1].boundary.plot(ax=axi, color="white", linewidth=0.4)
        axi.set_xlim(x0, x1)
        axi.set_ylim(y0, y1)
        axi.set_axis_on()
        axi.set_xticks([])
        axi.set_yticks([])
        axi.grid(False)
        for s in axi.spines.values():
            s.set_visible(True)
            s.set_color(bc.INK)
        axi.set_title("London", fontsize=12, loc="left", pad=3)
    if colours is not None:
        handles = []
        for k in colours:
            lab = labels.get(k, str(k)) if labels else str(k)
            if k == HATCH_VALUE:
                handles.append(Patch(facecolor="white", edgecolor=bc.INK, hatch="////", label=lab))
            else:
                handles.append(Patch(facecolor=colours[k], label=lab))
        if ns_mask is not None:
            handles.append(Patch(facecolor=bc.NS_GREY, label="not significant"))
        fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.01, 0.035 * 7.0 / fig_h), frameon=False,
                   fontsize=11.5, ncol=3, handlelength=1.2, columnspacing=0.8, handletextpad=0.4)
    else:
        cax = fig.add_axes([0.06, 0.145 * 7.0 / fig_h, 0.55, 0.022 * 7.0 / fig_h])
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(cbar_label or column, fontsize=12)
        if cbar_ticks is not None:
            cb.set_ticks(cbar_ticks)
        cb.outline.set_visible(False)
        if ns_mask is not None:
            fig.patches.append(Rectangle((0.66, 0.145 * 7.0 / fig_h), 0.04, 0.022 * 7.0 / fig_h,
                                         transform=fig.transFigure, facecolor=bc.NS_GREY, edgecolor="none"))
            fig.text(0.72, 0.156 * 7.0 / fig_h, "not significant\n(corrected t)", fontsize=10.5,
                     color=bc.INK2, va="center")
    if note:
        fig.text(0.01, 0.004, note, fontsize=10, color=bc.INK2, va="bottom", wrap=True)
    return bc.save_fig(fig, name)


def symmetric_norm(values, q=0.98):
    import matplotlib as mpl

    v = np.nanquantile(np.abs(values), q)
    return mpl.colors.TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v)


def repel(ys, min_gap):
    """Nudge label y-positions apart (keeps order) for direct line labels."""
    ys = np.asarray(ys, dtype=float)
    order = np.argsort(ys)
    out = ys.copy()
    for _ in range(200):
        moved = False
        s = out[order]
        for k in range(1, len(s)):
            if s[k] - s[k - 1] < min_gap:
                d = (min_gap - (s[k] - s[k - 1])) / 2
                s[k - 1] -= d
                s[k] += d
                moved = True
        out[order] = s
        if not moved:
            break
    return out
