# FEM Topologie Optimierer

**Abschlussprojekt  Software Design WS 2025/26**

**Autoren:** Urs Müller & Jan Vogel

**Betreuung:** Matthias Panny

**Institution:** MCI Innsbruck

---

## Projektübersicht

Dieses Projekt ist eine **Web-Anwendung von Streamlit** zur **Topologieoptimierung von 2D- und 3D-Strukturen** mithilfe der **Finite-Elemente-Methode**.

Ziel ist es, innerhalb eines vordefinierten Bauraums optimale Lastpfade zu berechnen und ineffizientes Material iterativ zu entfernen. Das Ergebnis ist eine **gewichtsoptimierte Struktur**, die sich direkt **visualisieren** und **für den 3D Druck exportieren** lässt.

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

![User Interface](images/UI.png)


### 5. Bedienung der UI

Die Bedienung der Webseite ist prozessorientiert strukturiert und folgt dem Arbeitsablauf:

1. **Definition der Geometrie (Tab „Geometrie“):** Initial erfolgt die Auswahl der Dimensionalität (2D oder 3D). Anschließend werden die globalen Abmessungen (Breite, Höhe, Tiefe), die zugrundeliegenden Materialeigenschaften (E-Modul) sowie die gewünschte Netzauflösung festgelegt. Durch die Bestätigung wird der Bauraum angegeben und das Gitter generiert.

2. **Festlegung der Randbedingungen (Tab „Randbedingungen“):** In diesem Schritt wird das statische System definiert. Kinematische Randbedingungen (Fest- und Loslager) sowie eingeprägte Lasten können mittels Koordinaten im System platziert werden. Alternativ steht eine Auswahl vordefinierter Standardlastfälle (wie z.B. der MBB-Balken oder eine 3D-Brückenstruktur) zur Verfügung, um Evaluierungsszenarien zeiteffizient zu laden.

3. **Parametrisierung und Optimierung (Tab „Optimierung“):** Hier erfolgt die Konfiguration des Optimierungsalgorithmus. Die steuernden Parameter umfassen die anvisierte Zielmasse (relativ zum Ausgangsvolumen), den Faktor der Verformungsdarstellung, sowie die Schrittweite des Materialabtrags pro Iterationsschritt. Durch das Auslösen des Startbefehls wird der FEM-Solver initiiert und der iterative Berechnungsprozess zur Topologieoptimierung ausgeführt. Der Fortschritt kann live beobachtet werden.

4. **Ergebnisanalyse und Export:** Nach erfolgreichem Abschluss der Berechnung stellt die Anwendung verschiedene Exportmöglichkeiten zur Verfügung. Die Topologie kann zur weiteren Dokumentation als Grafik (.png), als Animation (.gif) oder als Volumenmodell (.stl) für den 3D-Druck oder zur weiteren bearbeitung exportiert werden.

---

##  Minimalanforderungen 

Die Kernanforderungen der Aufgabenstellung wurden vollständig **objektorientiert in Python** umgesetzt:

* **Parametrische Geometrie**
  Definition von Breite, Höhe, Tiefe und Rasterauflösung (Massepunkte & Federn)

* **Randbedingungen**
  Festlager, Loslager und gerichtete Punktlasten

* **FEM-Solver**
  Aufbau der globalen Steifigkeitsmatrix und Lösung des linearen Gleichungssystems:
  $$K \cdot \vec{u} = \vec{F}$$

* **Optimierungsschleife**
  Berechnung der Verformungsenergie und iteratives Entfernen ineffizienter Knoten bis zur Zielmasse

* **Visualisierung**
  Darstellung von Geometrie, Randbedingungen und optimiertem Ergebnis

  ![2D Ansicht](images/2D_MBB-Balken.png)

* **Modell Speicher (JSON)**

  Die Struktur inkl. Randbedingungen kann jederzeit als JSON-Datei (`storage.json`) gespeichert und geladen werden.

---

## Erweiterungen 

### 1. 2D & 3D Strukturen

* Gittergenerierung als 2D-Fläche oder echtes 3D-Volumen
* Vollständig erweiterter Solver für Raumdiagonalen und 3D Transformationen

![3D Ansciht](images/3D_Ansicht.png)

### 2. Schnellerer Solver (SciPy Sparse)

* Verwendung von `scipy.sparse.csr_matrix`
* Lösung über `spsolve`
* Reduktion von Speicherbedarf und Rechenzeit

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

![Topologieoptimierung MBB-Balken](images/GIF_MBB-Balken.gif)

### 7. Symmetrie-Modus

* Implementierung eines zuschaltbaren Symmetrie-Modus.
* Der Materialabtrag wird über die Mittelachse gespiegelt, um symmetrische Bauteile zu erzeugen.

### 8. Klassische Vorlagen

* Vordefinierte Szenarien:

  * MBB-Balken
  * Kragarm
  * 3D-Brücke
  * Zentraler Druckstab

### 9. Heatmap-Visualisierung der Spannungen

* Über die Verformungsdarstellung hinaus werden die auftretenden Spannungen in den Elementen während und nach der Optimierung farblich als Heatmap visualisiert.

### 10. Reale Materialparameter

* Anstelle der vereinfachten Einheitssteifigkeiten aus der Aufgabenstellung k = 1 N/m erlaubt die Anwendung die Definition realer Werkstoffe.
* Integrierte Datenbank für Standardmaterialien (z. B. Baustahl S235, Aluminium, Holz) mit automatischer Umrechnung des materialspezifischen E-Moduls in die Elementsteifigkeitsmatrix.

### 11. Interaktives Lastfall-Management

* Alle gesetzten Randbedingungen und Lasten werden als übersichtliche Liste mit ihren exakten Koordinaten und Werten angezeigt und können einzeln oder gemeinsam gelöscht werden.

---

## Projektstruktur

```
├── User_Interface.py      # Streamlit-Webinterface
├── model.py               # FEM-Datenstrukturen 
├── solver.py              # Sparse FEM-Solver (SciPy)
├── UI_utils.py            # Visualisierung, Lasten, STL/GIF Export, Hilfsfunktionen
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
Die erzeugten STL Dateien werden analog dazu ausschließlich im **ASCII-STL-Format** generiert.

LLM´s wurden als Hilfsmittel für den Projektaufbau verwendet.
