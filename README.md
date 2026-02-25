# FEM Topologie Optimierer

**Abschlussprojekt  Software Design WS 2025/26**

**Autoren:** Urs Müller & Jan Vogel

**Betreuung:** Matthias Panny

**Institution:** MCI Innsbruck

---

## Projektübersicht

Dieses Projekt ist eine interaktive **Web-Anwendung auf Basis von Streamlit** zur **Topologieoptimierung von 2D- und 3D-Strukturen** mithilfe der **Finite-Elemente-Methode**.

Ziel ist es, innerhalb eines vordefinierten Bauraums optimale Lastpfade zu berechnen und ineffizientes Material iterativ zu entfernen. Das Ergebnis ist eine **gewichtsoptimierte Struktur**, die sich direkt **visualisieren** und **für den 3D-Druck exportieren** lässt.

---

##  Installation & Start

### 1. Repository klonen

```bash
git clone https://github.com/JanVogel004/Abschlussprojekt_Mueller_Vogel
cd Abschlussprojekt_Mueller_Vogel
```

### 2. Virtuelle Umgebung erstellen

```bash
python -m venv venv
venv\Scripts\activate
```

### 3. Abhängigkeiten installieren

```bash
pip install -r requirements.txt
```

### 4. Anwendung starten

```bash
streamlit run User_Interface.py
```

Die Anwendung öffnet sich automatisch im Standard Webbrowser.

---

##  Minimalanforderungen 

Die Kernanforderungen der Aufgabenstellung wurden vollständig **objektorientiert in Python** umgesetzt:

* **Parametrische Geometrie**
  Definition von Breite, Höhe, Tiefe und Rasterauflösung (Massepunkte & Federn)

* **Randbedingungen**
  Festlager, Loslager und gerichtete Punktlasten

* **FEM-Solver**
  Aufbau der globalen Steifigkeitsmatrix und Lösung von
  ( K * U = F )

* **Optimierungsschleife**
  Berechnung der Verformungsenergie und iteratives Entfernen ineffizienter Elemente bis zur Zielmasse

* **Visualisierung**
  Darstellung von Geometrie, Randbedingungen und optimiertem Ergebnis

---

## Erweiterungen 

### 1. 2D & 3D Strukturen

* Gittergenerierung als 2D-Fläche oder echtes 3D-Volumen
* Vollständig erweiterter Solver für Raumdiagonalen und 3D Transformationen

### 2. Schnellerer Solver (SciPy Sparse)

* Verwendung von `scipy.sparse.csr_matrix`
* Lösung über `spsolve`
* Massive Reduktion von Speicherbedarf und Rechenzeit

### 3. STL-Export

* **2D-Modelle:** 2.5D-Extrusion zu druckbaren Balken
* **3D-Modelle:** Echte Volumenmodelle durch lokale Orthogonalbasen

### 4. Strecken- & Flächenlasten

* Implementierung von Streckenlasten
* Gleichmäßige Lastverteilung zur Vermeidung von Singularitäten

### 5. Interaktive 3D Visualisierung

* Plotly basierte Darstellung
* Zoomen, Drehen, Anzeigen von Kraftvektoren (`go.Cone`)

### 6. Animations Export (.gif)

* Speicherung des gesamten Optimierungsverlaufs
* Export als GIF zur Visualisierung des Materialabbaus

### 7. Modell Speicher (JSON)

* Persistente Speicherung kompletter Lastfälle
* Laden & Speichern über `storage.json`

### 8. Klassische Vorlagen (Templates)

* Vordefinierte Szenarien:

  * MBB-Balken
  * Kragarm
  * 3D-Brücke
  * Zentraler Druckstab

---

## Projektstruktur

```
├── User_Interface.py      # Streamlit-Webinterface
├── model.py               # FEM-Datenstrukturen 
├── solver.py              # Sparse FEM-Solver (SciPy)
├── UI_utils.py            # Visualisierung, Lasten, STL/GIF Export
├── model_storage.py       # JSON Speicherlogik
├── storage.json           # Gespeicherte Lastfälle
├── requirements.txt       # Python Abhängigkeiten
```


## Verwendete Hilfsmittel & Inspirationen

Für die Entwicklung dieses Projekts haben wir uns an mehreren bestehenden Open Source Projekten orientiert und unterschiedliche Hilfsmittel eingesetzt.

Als zentrale Inspiration für die Topologieoptimierung diente unter anderem das GitHub Repository  
https://github.com/zfergus/topopt, welches ebenfalls auf `scipy.sparse` für die effiziente Lösung großer linearer Gleichungssysteme setzt.  
Uns ist bewusst, dass der Umstieg von NumPy auf SciPy Sparse bei den in diesem Projekt verwendeten Modellgrößen noch keinen signifikanten Performancegewinn bringt. Der Einsatz erfolgte dennoch bewusst aus technischem Interesse und zur Vorbereitung auf größere Problemstellungen.

Zusätzlich wurden weitere GitHub Projekte analysiert, um unterschiedliche Lösungsansätze, Datenstrukturen und Implementierungsdetails zu vergleichen.

Für den **STL-Export** wurde konzeptionelle Inspiration aus folgendem Repository entnommen:  
https://github.com/thearn/stl_tools  
Die erzeugten STL Dateien werden analog dazu ausschließlich im **offiziellen ASCII-STL-Format** generiert.

Zur Fehlersuche und zum Debugging komplexer Logik z. B. bei FEM Assemblierung, 3D Transformationen und Exportfunktionen wurden **Large Language Models** eingesetzt.