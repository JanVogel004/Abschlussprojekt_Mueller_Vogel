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


# User Interface--------------------------------------------------------------------

st.set_page_config(page_title="FEM Optimierung", layout="wide")
st.title("FEM Optimierer")

st.subheader("Gespeicherte Modelle")

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

    col1, col2 = st.columns([1,3])

    with col1:
        if st.button("Modell laden"):

            params = models[selected_index]

            # Parameter in session_state laden
            st.session_state.dims = (params["width"], params["height"])
            st.session_state.e_modul = params["e_mod"]
            st.session_state.constraints = params["constraints"]
            st.session_state.use_symmetry = params["symmetry"]
            st.session_state.res = params["width"] / (params["nx"] - 1)
            st.session_state.name = params["name"]

            nx = params["nx"]
            nz = params["nz"]
            e_mod = params["e_mod"]

            # Struktur neu erzeugen

            s = Structure(dim=2)
            # Grundstruktur wieder erzeugen
            s.generate_rect_mesh(params["nx"], params["nz"], params["width"], params["height"])

            for f in s.federn: f.k = params["e_mod"] / st.session_state.res

            for m in s.massepunkte:
                m.active = True
                m.displacement[:] = 0.0
                m.force[:] = 0.0
                m.fixed[:] = False

            st.session_state.structure = s

            st.success("Modell geladen")
            st.rerun()

else:
    st.info("Noch keine gespeicherten Modelle vorhanden.")


with st.sidebar:
    st.header("1. Geometrie")
    # Startwerte für Breite, Höhe, Auflösung und Material
    width = st.number_input("Gesamtbreite (m)", value=40.0, step=1.0)
    height = st.number_input("Gesamthöhe (m)", value=15.0, step=1.0)
    res = st.slider("Auflösung (m)", 0.5, 5.0, 1.5)
    st.divider()
    mat_type = st.selectbox("Material", ["Baustahl S235", "Aluminium", "Holz", "Custom"])
    e_mod_map = {"Baustahl S235": 210000.0, "Aluminium": 70000.0, "Holz": 10000.0, "Custom": 1000.0}
    e_modul = st.number_input("E-Modul (N/mm²)", value=e_mod_map[mat_type])
    
    if st.button("Gitter neu erzeugen", type="primary"):
        nx = int(width / res); nx = nx + 1 if nx % 2 == 0 else nx   # Kontenanzahl, immer ungerade für Symmetrie
        nz = int(height / res); nz = nz + 1 if nz % 2 == 0 else nz
        s = Structure(dim=2)
        s.generate_rect_mesh(nx, nz, width, height)   # Gitter erzeugen basierend auf Breite, Höhe und Auflösung
        for f in s.federn: f.k = e_modul / res
        st.session_state.structure = s
        st.session_state.dims = (width, height)
        st.session_state.e_modul = e_modul # Model speichern
        st.session_state.last_result_fig = None
        st.session_state.res = res
        st.session_state.constraints = [] # Alle Randbedingungen zurücksetzen
        st.session_state.use_symmetry = False
        st.session_state.name = None
        st.rerun()


if st.session_state.structure:
    s, (w, h), e_mod = st.session_state.structure, st.session_state.dims, st.session_state.e_modul
    apply_constraints(s)

    st.subheader("Ausgangslage")
    
    # Plot der Ausgangslage mit Randbedingungen, Kräften und eventuell Symmetrielinie
    st.pyplot(plot_with_stresses(s, "Vorschau", e_mod, w, h, vis_factor=0, is_setup_view=True, draw_sym_line=st.session_state.use_symmetry))

    #tabs für Randbedingungen und Optimierung
    tab1, tab2 = st.tabs(["Randbedingungen", "Optimierung"])
    
    with tab1:
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Lager**")
            lx = st.number_input("Position X", 0.0, w, 0.0, key="lx_in", step=1.0)
            lz = st.number_input("Position Z", 0.0, h, h, key="lz_in", step=1.0)
            ltype = st.radio("Typ", ["Festlager (∆)", "Loslager (○)"], key="ltype_in")
            if st.button("Lager anwenden"):
                clean_type = "Festlager" if "Fest" in ltype else "Loslager"
                st.session_state.constraints.append({"type": clean_type, "x": lx, "z": lz})
                st.rerun()
        with c2:
            # Kraft mit Richtung und Betrag ändern // steps in 1.0 schritten
            st.markdown("**Kraft**")
            kx = st.number_input("Position X", 0.0, w, w/2, key="kx_in", step=1.0)
            kz = st.number_input("Position Z", 0.0, h, 0.0, key="kz_in", step=1.0)
            fv = st.number_input("Betrag (N)", value=5000.0, key="fv_in", step=500.0)
            fw = st.number_input("Winkel (°)", 0.0, 360.0, 270.0, step=45.0, key="fw_in")
            if st.button("Kraft anwenden"): 
                st.session_state.constraints.append({"type": "Kraft", "x": kx, "z": kz, "val": fv, "angle": fw}) 
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

    with tab2:
        # Checkbox Symetrie-Modus
        st.session_state.use_symmetry = st.checkbox("Symmetrie-Modus nutzen (Spiegelt Materialabtrag)", value=st.session_state.use_symmetry)

        target = st.slider("Ziel-Masse (%)", 5, 100, 70) / 100.0
        step = st.slider("Schrittweite", 0.01, 0.1, 0.02)
        vis = st.slider("Verformungs-Faktor", 1.0, 10.0, 1.0)
        
        st.session_state.name = st.text_input("Name des Modells", placeholder="z.B. Holzbalken mit mittiger Last", value = st.session_state.name if st.session_state.name else "")

        if st.button("Optimierung starten", type="primary", disabled=not st.session_state.name.strip()):
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
                fig = plot_with_stresses(c_struct, "Optimierung", e_mod, w, h, vis, 
                                    draw_sym_line=st.session_state.use_symmetry, 
                                    current_mass_pct=curr_m*100)
                plot_spot.pyplot(fig); plt.close(fig)

            st.session_state.last_result_fig = fig

            # Ergebnis kann als PNG heruntergeladen werden
            if st.session_state.last_result_fig:

                png_bytes = fig_to_png_bytes(st.session_state.last_result_fig)

                st.download_button(
                    label="Geometrie als PNG herunterladen",
                    data=png_bytes,
                    file_name="optimierte_geometrie.png",
                    mime="image/png"
                )