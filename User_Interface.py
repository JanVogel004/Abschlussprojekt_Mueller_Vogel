import streamlit as st
import numpy as np
import matplotlib
matplotlib.use('Agg') # matplotlib im Hintergrundmodus für Streamlit
import matplotlib.pyplot as plt
import io
from model import Structure
from model_storage import save_model, load_all_models, get_model_names
from UI_utils import *

# speichert bei neu laden den alten Stand
if 'structure' not in st.session_state:
    st.session_state.structure = None
if 'last_result_fig' not in st.session_state:
    st.session_state.last_result_fig = None
if 'use_symmetry' not in st.session_state:
    st.session_state.use_symmetry = False
if 'constraints' not in st.session_state:
    st.session_state.constraints = []
if 'res' not in st.session_state:
    st.session_state.res = None
if 'name' not in st.session_state:
    st.session_state.name = None
if 'target_mass' not in st.session_state:
    st.session_state.target_mass = 0.60
if 'step_size' not in st.session_state:
    st.session_state.step_size = 0.04


# User Interface--------------------------------------------------------------------

st.set_page_config(page_title="FEM Optimierung", layout="wide")
st.title("FEM Optimierer")

# --- CSS HACK: "Press Enter to apply" verstecken ---
hide_streamlit_instructions = """
<style>
    /* Versteckt die kleinen Hilfetexte unter den Eingabefeldern */
    div[data-testid="InputInstructions"] {
        display: none !important;
    }
</style>
"""
st.markdown(hide_streamlit_instructions, unsafe_allow_html=True)
#--------------------------------------------------------------


with st.sidebar:
    st.header("Gespeicherte Modelle")

    models = load_all_models()

    if models:

        # Anzeige-Text für jedes Modell erzeugen
        selected_index = st.selectbox(
        "Modell auswählen",
        range(len(models)),
        format_func=lambda i: f" {models[i]['name']} ({models[i]['date']})"
        )
        m = models[selected_index]
        # Details zum Modell anzeigen
        st.markdown(
            f"""
        **Geometrie:** {m['width']} × {m['height']} m  
        **Mesh:** {m['nx']} × {m['nz']} Knoten  
        **Material:** E = {m['e_mod']} N/mm²  
        **Randbedingungen:** {len(m['constraints'])} Einträge  
        """
        )
    else:
        st.info("Noch keine gespeicherten Modelle vorhanden.")

    if st.button("Modell laden"):

        params = models[selected_index]

        # 1. Parameter in session_state laden
        is_3d = "depth" in params and params["depth"] > 0
        
        if is_3d:
            st.session_state.dims = (params["width"], params["height"], params["depth"])
        else:
            st.session_state.dims = (params["width"], params["height"])

        st.session_state.e_modul = params["e_mod"]
        st.session_state.constraints = params["constraints"]
        st.session_state.use_symmetry = params["symmetry"]
        st.session_state.res = params["res"]
        st.session_state.name = params["name"]
        st.session_state.target_mass = params.get("target_mass", 0.60)
        st.session_state.step_size = params.get("step", 0.04)
        st.session_state.nx = params["nx"]
        st.session_state.nz = params["nz"]
        st.session_state.ny = params.get("ny", 1)

        # 2. Struktur erzeugen
        if is_3d:
            s = Structure(dim=3)
            s.generate_box_mesh(params["nx"], params["ny"], params["nz"], 
                                params["width"], params["depth"], params["height"])
        else:
            s = Structure(dim=2)
            s.generate_rect_mesh(params["nx"], params["nz"], params["width"], params["height"])

        # 3. Materialsteifigkeit und Status zurücksetzen
        for f in s.federn: f.k = params["e_mod"] / st.session_state.res
        
        for m in s.massepunkte:
            m.active = True
            m.displacement[:] = 0.0
            m.force[:] = 0.0
            m.fixed[:] = False

        st.session_state.structure = s
        st.session_state.last_result_fig = None # Altes Bild löschen
        st.success("Modell geladen!")
        st.rerun()



if st.session_state.structure:
    s, dims, e_mod = st.session_state.structure, st.session_state.dims, st.session_state.e_modul
    w, h = dims[0], dims[1]
    d = dims[2] if len(dims) > 2 else 0.0  # Holt die Tiefe, wenn 3D
    apply_constraints(s)

    st.subheader("Ausgangslage")
    if s.dim == 2:
        st.pyplot(plot_with_stresses(s, "Vorschau", e_mod, w, h, vis_factor=0, is_setup_view=True, draw_sym_line=st.session_state.use_symmetry))
    else:
        fig_3d = plot_3d_structure(s, "Vorschau", e_mod, vis_factor=0, is_setup_view=True)
        st.plotly_chart(fig_3d, use_container_width=True)

else:
    st.info("Willkommen! Bitte erzeuge im Tab 'Geometrie' zuerst ein Gitter.")

st.divider()

#tabs für Randbedingungen und Optimierung
tab1, tab2, tab3 = st.tabs(["Geometrie", "Randbedingungen", "Optimierung"])


with tab1:
    dim_mode = st.radio("Dimension", ["2D", "3D"], horizontal=True)
    with st.form("geometrie_form"):
        # Startwerte für Breite, Höhe, Tiefe, Auflösung und Material
        col_w, col_h, col_d = st.columns(3)
        with col_w:
            width = st.number_input("Breite X (m)", value=40.0, step=1.0)
        with col_h:
            height = st.number_input("Höhe Z (m)", value=15.0, step=1.0)
        with col_d:
            # Tiefe ist nur in 3D anwählbar
            depth = st.number_input("Tiefe Y (m)", value=10.0, step=1.0, disabled=(dim_mode == "2D"))
        res = st.slider("Auflösung (m)", 0.1, 4.0, 1.5, step=0.1)
        st.divider()

        mat_type = st.selectbox("Material", ["Baustahl S235", "Aluminium", "Holz", "Custom"])
        e_mod_map = {"Baustahl S235": 210000.0, "Aluminium": 70000.0, "Holz": 10000.0, "Custom": 1000.0}
        e_modul = st.number_input("E-Modul (N/mm²)", value=e_mod_map[mat_type])

        if st.form_submit_button("Gitter neu erzeugen", type="primary"):
            nx = int(width / res); nx = nx + 1 if nx % 2 == 0 else nx   # Kontenanzahl, immer ungerade für Symmetrie
            nz = int(height / res); nz = nz + 1 if nz % 2 == 0 else nz
            # Weiche für 2D / 3D
            if dim_mode == "2D":
                s = Structure(dim=2)
                s.generate_rect_mesh(nx, nz, width, height)
                st.session_state.dims = (width, height) # 2 Werte speichern
            else:
                ny = int(depth / res); ny = ny + 1 if ny % 2 == 0 else ny
                s = Structure(dim=3)
                s.generate_box_mesh(nx, ny, nz, width, depth, height)
                st.session_state.dims = (width, height, depth) # 3 Werte speichern
                st.session_state.ny = ny
            
            st.session_state.nx = nx
            st.session_state.nz = nz

            for f in s.federn: f.k = e_modul / res
            st.session_state.structure = s

            st.session_state.e_modul = e_modul # Model speichern
            st.session_state.last_result_fig = None
            st.session_state.res = res
            st.session_state.constraints = [] # Alle Randbedingungen zurücksetzen
            st.session_state.use_symmetry = False
            st.session_state.name = None
            st.rerun()

with tab2:
    if not st.session_state.structure:
        st.warning("Bitte erstelle zuerst ein Gitter im Tab 'Geometrie'!")
    else:
        s, dims, e_mod = st.session_state.structure, st.session_state.dims, st.session_state.e_modul
        w, h = dims[0], dims[1]
        d = dims[2] if len(dims) > 2 else 0.0  # Holt die Tiefe, wenn 3D

        c1, c2 = st.columns(2)
        
        with c1:
            st.markdown("**Lager**")
            with st.form("lager_form"):
                lx = st.number_input("Position X", 0.0, w, 0.0, key="lx_in", step=1.0)
                lz = st.number_input("Position Z", 0.0, h, h, key="lz_in", step=1.0)
                ly = st.number_input("Position Y (Tiefe)", 0.0, d, 0.0, step=1.0, disabled=(s.dim==2))
                ltype = st.radio("Typ", ["Festlager (∆)", "Loslager (○)"], key="ltype_in")
                
                if st.form_submit_button("Lager anwenden"):
                    clean_type = "Festlager" if "Fest" in ltype else "Loslager"
                    st.session_state.constraints.append({"type": clean_type, "x": lx, "z": lz, "y": ly})
                    st.rerun()
        with c2:
            # Kraft mit Richtung und Betrag ändern // steps in 1.0 schritten
            st.markdown("**Kraft**")
            with st.form("kraft_form"):
                kx = st.number_input("Position X", 0.0, w, w/2, key="kx_in", step=1.0)
                kz = st.number_input("Position Z", 0.0, h, 0.0, key="kz_in", step=1.0)
                ky = st.number_input("Position Y (Tiefe)", 0.0, d, d/2, step=1.0, disabled=(s.dim==2))
                fv = st.number_input("Betrag (N)", value=5000.0, key="fv_in", step=500.0)
                if s.dim == 2:
                    fw = st.number_input("Winkel (°)", 0.0, 360.0, 270.0, step=45.0)
                else:
                    st.info("In 3D wirkt die Kraft standardmäßig senkrecht nach unten.")
                
                if st.form_submit_button("Kraft anwenden"): 
                    constraint = {"type": "Kraft", "x": kx, "z": kz, "y": ky, "val": fv}
                    if s.dim == 2: constraint["angle"] = fw
                    st.session_state.constraints.append(constraint) 
                    st.rerun()

        st.divider()
        # Liste der Randbedingungen
        st.write("**Aktuelle Lasten & Lager:**")

        if st.session_state.constraints:
            for i, c in enumerate(st.session_state.constraints):
                col_text, col_del = st.columns([4, 1])
                
                if c['type'] == "Kraft":
                    # Anzeige für KRAFT
                    info_text = f"**Kraft**: {c['val']} N ({c.get('angle', 270)}°) an Position ({c['x']}, {c['z']})"
                else:
                    # Anzeige für LAGER
                    symbol = "∆" if "Fest" in c['type'] else "○"
                    info_text = f"**{c['type']}** {symbol} an Position ({c['x']}, {c['z']})"

                col_text.markdown(info_text)
                
                # Löschen-Button
                if col_del.button("Löschen", key=f"delete_{i}"):
                    st.session_state.constraints.pop(i) # Eintrag entfernen
                    st.rerun() # Seite neu laden
        else:
            st.info("Die Liste ist leer. Füge oben Lager oder Kräfte hinzu.")
        st.divider()

        col1, col2 = st.columns([1,1])

        # Reset Button
        with col1:
            if st.button("Alles löschen (Reset)"):
                st.session_state.constraints = [] # Liste leeren!
                apply_constraints(s) # Modell leeren
                st.rerun()
        # Zu Optimierung wechseln
        with col2:
            st.text("Nun weiter zur Optimierung!")

with tab3:
    if not st.session_state.structure:
        st.warning("Bitte erstelle zuerst ein Gitter im Tab 'Geometrie'!")
    else:
        s, dims, e_mod = st.session_state.structure, st.session_state.dims, st.session_state.e_modul
        w, h = dims[0], dims[1]
        d = dims[2] if len(dims) > 2 else 0.0  # Holt die Tiefe, wenn 3D

        with st.form("optimierung_form"):
            # Checkbox Symetrie-Modus
            st.session_state.use_symmetry = st.checkbox("Symmetrie-Modus nutzen (Spiegelt Materialabtrag)", value=st.session_state.use_symmetry)

            default_target = int(st.session_state.target_mass * 100)
            target = st.slider("Ziel-Masse (%)", 5, 100, default_target) / 100.0
            step = st.slider("Schrittweite", 0.01, 0.1, st.session_state.step_size)
            vis = st.slider("Verformungs-Faktor", 1.0, 10.0, 1.0)
            
            # Namensfeld
            name_input = st.text_input("Name des Modells", value=st.session_state.name if st.session_state.name else "")

            col_btn1, col_btn2 = st.columns(2)
            with col_btn1:
                save_clicked = st.form_submit_button("Modell speichern")
            with col_btn2:
                start_clicked = st.form_submit_button("Optimierung starten", type="primary")

        if save_clicked or start_clicked:
            # Daten im Session State aktualisieren
            st.session_state.name = name_input.strip()
            st.session_state.target_mass = target
            st.session_state.step_size = step

            if save_clicked:
                clean_name = st.session_state.name.strip()
                
                if not clean_name:
                    st.error("Bitte gib einen Namen ein, um zu speichern!")
                elif clean_name in get_model_names():
                    st.error(f"Der Name '{clean_name}' existiert bereits! Bitte wähle einen anderen.")
                else:
                    save_current_model(target, step, clean_name)
                    st.success(f"Modell '{clean_name}' wurde gespeichert.")
                    st.rerun()    
            
            if start_clicked:    
                apply_constraints(s)
                
                # Testen ob Struktur stabil ist, bevor optimiert wird
                stable, message = s.stable_test()
                if not stable:
                    st.error(f"Die Struktur ist instabil: {message} :(")
                    st.stop()
                st.success("Die Struktur ist stabil! Starte Optimierung...")

                # Parameter speichern für später, wenn es ein Neues Modell ist
                if (st.session_state.name.strip() not in get_model_names()):
                    save_current_model(target, step, st.session_state.name.strip())

                # Volle Struktur nutzen
                c_struct = s 
                for m in c_struct.massepunkte: m.active, m.displacement[:] = True, 0
                
                # Symmetrie-Partner vorberechnen
                partner_map = {} 
                if st.session_state.use_symmetry:
                    nodes = c_struct.massepunkte
                    for i, n in enumerate(nodes):
                        if n.coords[0] < w/2 - 0.01: # Links
                            target_x = w - n.coords[0]
                            target_z = n.coords[1]
                            dists = np.linalg.norm([p.coords - np.array([target_x, target_z]) for p in nodes], axis=1)
                            best_match_idx = np.argmin(dists)
                            if dists[best_match_idx] < 0.1: 
                                partner_map[i] = best_match_idx

                plot_spot = st.empty()
                curr_m = 1.0
                
                while curr_m > target:  # Rechne bis Zielmasse erreicht ist
                    if c_struct.solve() is None: break  # Verschiebung berechnen
                    c_struct.calculate_strain_energy()  # Energie berechnen
                    
                    # 1. Energie mitteln (links und rechts) für Symetrie
                    if st.session_state.use_symmetry:
                        for idx_left, idx_right in partner_map.items():
                            avg_energy = (c_struct.massepunkte[idx_left].strain_energy + c_struct.massepunkte[idx_right].strain_energy) / 2
                            c_struct.massepunkte[idx_left].strain_energy = avg_energy
                            c_struct.massepunkte[idx_right].strain_energy = avg_energy
                    
                    # 2. Schutz für Massepunkte mit Kraft oder Randbedingung: Sehr hohe Energie zuweisen, damit sie nicht gelöscht werden
                    for m in c_struct.massepunkte:
                        if np.linalg.norm(m.force) > 0 or np.any(m.fixed):
                            m.strain_energy = 1e15
                    
                    # 3. Löschen von Massepunkten mit geringster Energie je nach Schrittweite 1 bis 10% der Massepunkte
                    c_struct.remove_inefficient_nodes(max(target, curr_m - step))

                    # 4. Konstrolle der Symetrie
                    if st.session_state.use_symmetry:
                        for idx_left, idx_right in partner_map.items():
                            # Wenn einer der beiden inaktiv ist, muss der andere es auch sein
                            if not c_struct.massepunkte[idx_left].active or not c_struct.massepunkte[idx_right].active:
                                c_struct.massepunkte[idx_left].active = False
                                c_struct.massepunkte[idx_right].active = False

                    # 5. Plotten der aktuellen Struktur mit Spannungen
                    curr_m = len([m for m in c_struct.massepunkte if m.active]) / len(c_struct.massepunkte)
                    
                    if s.dim == 2:
                        fig = plot_with_stresses(c_struct, "Optimierung", e_mod, w, h, vis, 
                                            draw_sym_line=st.session_state.use_symmetry, 
                                            current_mass_pct=curr_m*100)
                        plot_spot.pyplot(fig); plt.close(fig)
                        st.session_state.last_result_fig = fig
                    else:
                        fig = plot_3d_structure(c_struct, "Optimierung", e_mod, vis, 
                                            current_mass_pct=curr_m*100)
                        plot_spot.plotly_chart(fig, use_container_width=True)
                        st.session_state.last_result_fig = fig

                st.session_state.last_result_fig = fig

    # Ergebnis kann als PNG heruntergeladen werden
    if st.session_state.last_result_fig:
        if s.dim == 3:
            # Manueller Export für Plotly
            btn_data = plotly_to_png_bytes(st.session_state.last_result_fig)
        else:
            btn_data = fig_to_png_bytes(st.session_state.last_result_fig)

        st.download_button(
            label="3D-Ansicht als PNG speichern",
            data=btn_data,
            file_name="3d_optimierung.png",
            mime="image/png"
        )