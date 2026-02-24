import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import io
from model_storage import save_model  # Wird für save_current_model gebraucht

#Hilfsfunktionen

# Funktion um Parameter eines Modells zu speichern
def save_current_model(target, step, name):

    width, height = st.session_state.dims

    save_model(
        name=name,

        width=width,
        height=height,

        nx=int(width / st.session_state.res),
        nz=int(height / st.session_state.res),
        res=st.session_state.res,

        e_mod=st.session_state.e_modul,

        constraints=st.session_state.constraints,

        symmetry=st.session_state.use_symmetry,

        target_mass=target,
        step=step
    )

# Funktion um Figur herunterzuladen
def fig_to_png_bytes(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    buf.seek(0)
    return buf

# Plot Gitter Mash--------------------------------------------------------------------
def plot_with_stresses(structure, title, e_mod, w_orig, h_orig, vis_factor, is_setup_view=False, current_mass_pct=None, draw_sym_line=False):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlim(-w_orig * 0.05, w_orig * 1.05)
    ax.set_ylim(-h_orig * 0.05, h_orig * 1.05)
    
    full_title = f"{title} (Material: {current_mass_pct:.1f}%)" if current_mass_pct is not None else title

    # Symmetrielinie
    if draw_sym_line:
        ax.axvline(x=w_orig/2, color='#0000FF', linestyle='--', alpha=0.8, lw=2, zorder=1)

    # Spannungen der Federn berechnen
    all_sigmas = []
    if not is_setup_view:
        for f in structure.federn:
            if f.massepunkt_i.active and f.massepunkt_j.active:
                L0 = np.linalg.norm(f.massepunkt_j.coords - f.massepunkt_i.coords)
                if L0 > 1e-9:
                    p1 = f.massepunkt_i.coords + f.massepunkt_i.displacement[:2]
                    p2 = f.massepunkt_j.coords + f.massepunkt_j.displacement[:2]
                    L1 = np.linalg.norm(p2 - p1)
                    all_sigmas.append(abs((L1 - L0) / L0) * e_mod)
    
    max_s = np.percentile(all_sigmas, 95) if all_sigmas and len(all_sigmas) > 5 else 1.0
    if max_s < 1e-3: max_s = 1.0

    cmap = plt.get_cmap('plasma')

    # Mash erzeugen
    for f in structure.federn:
        if is_setup_view or (f.massepunkt_i.active and f.massepunkt_j.active):
            p1 = f.massepunkt_i.coords + f.massepunkt_i.displacement[:2] * vis_factor
            p2 = f.massepunkt_j.coords + f.massepunkt_j.displacement[:2] * vis_factor
            
            # Layer 1: Skelett
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color='black', lw=0.5, alpha=0.2, zorder=1)

            color = '#A0A0A0'
            lw = 0.3
            alpha = 0.5

            if not is_setup_view:
                L0 = np.linalg.norm(f.massepunkt_j.coords - f.massepunkt_i.coords)
                sigma = 0
                if L0 > 1e-9:
                    p1_real = f.massepunkt_i.coords + f.massepunkt_i.displacement[:2]
                    p2_real = f.massepunkt_j.coords + f.massepunkt_j.displacement[:2]
                    sigma = abs((np.linalg.norm(p2_real-p1_real) - L0)/L0 * e_mod)
                
                color = cmap(min(sigma / max_s, 1.0))
                lw = 0.5 + (min(sigma / max_s, 1.0) * 3.0)
                alpha = 0.9

            # Layer 2: Farbe
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, linewidth=lw, alpha=alpha, zorder=2)
    
    for m in structure.massepunkte:
        if is_setup_view or m.active:
            pos = m.coords + m.displacement[:2] * vis_factor
            if np.all(m.fixed): ax.plot(pos[0], pos[1], '^', color='red', ms=12, zorder=10)
            elif np.any(m.fixed): ax.plot(pos[0], pos[1], 'o', mfc='none', mec='orange', mew=2, ms=10, zorder=10)
            if np.linalg.norm(m.force) > 0:
                scale = (h_orig * 0.15) / np.linalg.norm(m.force)
                dx = m.force[0] * scale
                dy = m.force[1] * scale
                ax.arrow(pos[0], pos[1], dx, dy, head_width=w_orig*0.02, color='#00FF00', lw=2, zorder=20)
    
    if not is_setup_view:
        sm = plt.cm.ScalarMappable(cmap=plt.get_cmap('plasma'), norm=plt.Normalize(vmin=0, vmax=max_s))
        plt.colorbar(sm, ax=ax, label="Spannung [N/mm²]")

    ax.invert_yaxis(); ax.set_aspect('equal'); ax.set_title(full_title)
    return fig

def apply_constraints(struct):
    # 1. Alles zurücksetzen
    for m in struct.massepunkte: 
        m.fixed[:] = False
        m.force[:] = 0.0
    
    # 2. Liste aus Streamlit durchgehen
    for c in st.session_state.constraints:
        target = np.array([c['x'], c['z']]) # Sucht nächsten Knoten zum angegebenen Punkt
        node = min(struct.massepunkte, key=lambda m: np.linalg.norm(m.coords - target))
        
        if c['type'] == "Festlager":    # setzt Lager oder Kraft
            node.fixed[:] = True
        elif c['type'] == "Loslager": 
            node.fixed[1] = True 
        elif c['type'] == "Kraft":
            w_grad = c.get('angle', 270.0)
            w_rad = np.radians(w_grad)
            F = c['val']
            #mehrere Kräfte gleichzeitig möglich
            node.force[0] += F * np.cos(w_rad) 
            node.force[1] += -1.0 * F * np.sin(w_rad)

def create_stl_from_2D_structure(structure, thickness=1.0, beam_width=1.0):
    # STL Datei generieren
    # STL ist eine lange Liste von Facetten, jede definiert durch 3 Punkte und einen Normalenvektor
    stl_lines = ["solid FEM_Export"]
    
    def add_facet(v1, v2, v3):
        # Einfacher Normalenvektor
        n = np.cross(v2 - v1, v3 - v1)
        n = n / np.linalg.norm(n) if np.linalg.norm(n) > 0 else n
        stl_lines.append(f"  facet normal {n[0]} {n[1]} {n[2]}")
        stl_lines.append("    outer loop")
        for v in [v1, v2, v3]:
            stl_lines.append(f"      vertex {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}")
        stl_lines.append("    endloop")
        stl_lines.append("  endfacet")

    for f in structure.federn:
        if f.massepunkt_i.active and f.massepunkt_j.active:
            p1, p2 = f.massepunkt_i.coords, f.massepunkt_j.coords
            vec = p2 - p1
            L = np.linalg.norm(vec)
            if L < 1e-6: continue
            
            # Ausrichtung für Balkendicke
            dir_vec = vec / L
            perp = np.array([-dir_vec[1], dir_vec[0]]) * (beam_width / 2)
            z_off = thickness / 2
            
            # 8 Eckpunkte des 3D-Balkens 
            nodes = [
                np.array([p1[0]-perp[0], p1[1]-perp[1],  z_off]), # 0
                np.array([p1[0]+perp[0], p1[1]+perp[1],  z_off]), # 1
                np.array([p1[0]+perp[0], p1[1]+perp[1], -z_off]), # 2
                np.array([p1[0]-perp[0], p1[1]-perp[1], -z_off]), # 3
                np.array([p2[0]-perp[0], p2[1]-perp[1],  z_off]), # 4
                np.array([p2[0]+perp[0], p2[1]+perp[1],  z_off]), # 5
                np.array([p2[0]+perp[0], p2[1]+perp[1], -z_off]), # 6
                np.array([p2[0]-perp[0], p2[1]-perp[1], -z_off])  # 7
            ]
            # 12 Facetten für einen geschlossenen Quader
            faces = [[0,1,5],[0,5,4], [1,2,6],[1,6,5], [2,3,7],[2,7,6], [3,0,4],[3,4,7], [4,5,6],[4,6,7], [3,2,1],[3,1,0]]
            for face in faces:
                add_facet(nodes[face[0]], nodes[face[1]], nodes[face[2]])
                
    stl_lines.append("endsolid FEM_Export")
    return "\n".join(stl_lines).encode('utf-8')