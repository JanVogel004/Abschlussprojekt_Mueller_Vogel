import numpy as np
from dataclasses import dataclass, field
from solver import solve as linear_solve 
from scipy.sparse import csr_matrix

#scipy.sparse spart Speicherplatz und Rechenzeit

@dataclass
class Massepunkt:
    id: int             # eindeutige ID
    coords: np.ndarray  # Pos. [x, y, (z)]
    fixed: np.ndarray   # Abhängigkeiten (z.B. Festlager) [True, False]
    force: np.ndarray   # Kraftvektor [Fx, Fy, (Fz)]
    active: bool = True # False -> gelöscht
    mass: float = 1.0   # Standardmasse 1kg
    displacement: np.ndarray = field(default_factory=lambda: np.zeros(3))   # Verschiebungsvektor
    strain_energy: float = 0.0  #Verformungsenergie

class Feder:
    def __init__(self, massepunkt_i: Massepunkt, massepunkt_j: Massepunkt, k_val: float = 1.0):
        self.massepunkt_i = massepunkt_i    # Verbundene Massepunkte i und j
        self.massepunkt_j = massepunkt_j
        self.k = k_val                      # Federsteifigkeit (N/m)

    def get_element_stiffness(self, dim: int) -> np.ndarray:
        """
        Berechnet Steifigkeits- und Transformationsmatrix für 2D und 3D!
        """
        # Richtungsvektor berechnen 
        vec = self.massepunkt_j.coords - self.massepunkt_i.coords
        length = np.linalg.norm(vec)
        
        if length == 0:
            return np.zeros((2*dim, 2*dim)) # Vermeidung von Division durch Null -> gibt Matrix mit Nullen zurück
            
        e_n = vec / length # Einheitsvektor
        
        # Transformationsmatrix O = e_n * e_n
        O = np.outer(e_n, e_n)
        
        # Basis-Steifigkeit
        K_base = self.k * np.array([[1, -1], [-1, 1]])
        
        # Kronecker-Produkt für die globale Orientierung
        # Das erzeugt autom. die richtige Matrixgröße (4x4 für 2D, 6x6 für 3D)
        K_local = np.kron(K_base, O)
        
        return K_local  # Steifigkeits- mit Transformationsmatrix

class Structure:
    def __init__(self, dim: int = 2):
        self.massepunkte = []    # Alle Massepunkte der Struktur
        self.federn = []         # Alle Federn der Struktur
        self.dim = dim           # 2D oder 3D

    def add_Massepunkt(self, x, z, y=0.0):
        # Erstellt Massepunkt passend zur Dimension
        if self.dim == 2:
            # Im 2D-Fall speichern wir [x, z]
            coords = np.array([x, z])
        else:
            # Im 3D-Fall speichern wir [x, z, y]
            coords = np.array([x, z, y])
        # Initial keine Kraft, nicht fixiert
        fixed = np.zeros(self.dim, dtype=bool)
        force = np.zeros(self.dim)
        
        new_massepunkt = Massepunkt(len(self.massepunkte), coords, fixed, force)
        self.massepunkte.append(new_massepunkt)
        return new_massepunkt

    def add_Feder(self, idx_i, idx_j, k=1.0):
        # neue Feder zwischen Massepunkt i und j erstellen
        el = Feder(self.massepunkte[idx_i], self.massepunkte[idx_j], k)
        self.federn.append(el)
        
    def assemble_system(self):
            # Steifigkeitsmatrix K und Kraftvektor F aufbauen
            n_dof = len(self.massepunkte) * self.dim
            K_g = np.zeros((n_dof, n_dof))
            F_g = np.zeros(n_dof)
            
            # 1. Globale Steifigkeitsmatrix aufbauen
            for feder in self.federn:
                # Nur aktive Federn berücksichtigen!
                if not feder.massepunkt_i.active or not feder.massepunkt_j.active:
                    continue

                # Lokale Matrix holen (automatisch 2D oder 3D)
                k_local = feder.get_element_stiffness(self.dim)
                
                # IDs der beteiligten Knoten
                i = feder.massepunkt_i.id
                j = feder.massepunkt_j.id
                
                # Liste aller betroffenen Knotenindizes in der globalen Matrix:
                # Indizes für Knoten i (Start) und j (Ende)
                idx_i = list(range(i * self.dim, i * self.dim + self.dim))
                idx_j = list(range(j * self.dim, j * self.dim + self.dim))
                
                # Alle Indizes zusammen (das entspricht den Zeilen/Spalten in k_local)
                idxs = idx_i + idx_j
                
                # Wir addieren die lokale Matrix auf die Stellen in der globalen Matrix.
                K_g[np.ix_(idxs, idxs)] += k_local

            # 2. Kraftvektor aufbauen
            # Wir gehen alle Punkte durch und schauen, ob eine Kraft wirkt
            for node in self.massepunkte:
                # Wo startet der Eintrag für diesen Knoten im Vektor?
                start_idx = node.id * self.dim
                
                # Kraft eintragen (Fx, Fy, (Fz))
                F_g[start_idx : start_idx + self.dim] = node.force

            return K_g, F_g
    
    def assemble_stiffness_sparse(self):
        rows, cols, data = [], [], []
        n_dof = len(self.massepunkte) * self.dim
        # Nur aktive Federn berücksichtigen
        for feder in self.federn:
            if not feder.massepunkt_i.active or not feder.massepunkt_j.active:
                continue

            k_local = feder.get_element_stiffness(self.dim)

            i = feder.massepunkt_i.id
            j = feder.massepunkt_j.id
            # Alle Indizes der beteiligten Knoten in der globalen Matrix sammeln
            dofs = (
                list(range(i*self.dim, i*self.dim+self.dim)) +
                list(range(j*self.dim, j*self.dim+self.dim))
            )
            
            for a in range(2*self.dim):
                for b in range(2*self.dim):
                    val = k_local[a, b]
                    if abs(val) > 1e-12:
                        rows.append(dofs[a])
                        cols.append(dofs[b])
                        data.append(val)
        # Sparse-Matrix erstellen
        self.K_base = csr_matrix((data, (rows, cols)), shape=(n_dof, n_dof))

    def assemble_force_vector(self):
        # Kraftvektor F aufbauen
        F = np.zeros(len(self.massepunkte) * self.dim)
        for m in self.massepunkte:
            start = m.id * self.dim
            F[start:start+self.dim] = m.force
        return F
    
    def get_fixed_dofs(self):
        # Alle Freiheitsgrade sammeln, die fixiert oder inaktive Knoten sind
        fixed = []
        for m in self.massepunkte:
            for d in range(self.dim):
                idx = m.id * self.dim + d
                if m.fixed[d] or not m.active:
                    fixed.append(idx)
        return fixed
        
    def solve(self):
    # Steifigkeitsmatrix nur einmal aufbauen
        if not hasattr(self, "K_base"):
            self.assemble_stiffness_sparse()

        F = self.assemble_force_vector()
        fixed_dofs = self.get_fixed_dofs()

        u_vec = linear_solve(self.K_base, F, fixed_dofs)

        if u_vec is None:
            print("System ist singulär / nicht lösbar")
            return None

        for m in self.massepunkte:
            start = m.id * self.dim
            m.displacement = u_vec[start:start+self.dim]

        return u_vec

    def calculate_strain_energy(self):

        """
        Berechnet die Verformungsenergie für alle aktiven Federn
        und verteilt sie 50/50 auf die Massenpunkte.
        """

        # 1. Energie aller Knoten zurücksetzen
        for m in self.massepunkte:
            m.strain_energy = 0.0

        # 2. Über alle Federn iterieren
        for feder in self.federn:
            # Ignoriere inaktive Federn
            if not feder.massepunkt_i.active or not feder.massepunkt_j.active:
                continue

            # Vektoren und Matrizen holen
            k_local = feder.get_element_stiffness(self.dim)
            
            # Verschiebeungsvektor u für dieses Element bauen
            u_i = feder.massepunkt_i.displacement[:self.dim]
            u_j = feder.massepunkt_j.displacement[:self.dim]
            u_element = np.concatenate([u_i, u_j])

            # Matrix-Multiplikation: 0.5 * u^T * K * u
            # (u @ K) -> Vektor-Matrix-Multiplikation
            # dot(u) -> Skalarprodukt am Ende
            energy = 0.5 * np.dot(u_element @ k_local, u_element)

            # Energie 50/50 aufteilen
            feder.massepunkt_i.strain_energy += energy * 0.5
            feder.massepunkt_j.strain_energy += energy * 0.5

    def remove_inefficient_nodes(self, target_mass_percent: float):
        """
        Deaktiviert die Knoten mit dem geringsten Kraftfluss, bis die Zielmasse erreicht ist.
        target_mass_percent: z.B. 0.5 für 50% der Masse übrig lassen.
        """
        # Nur momentan aktive Knoten betrachten
        active_nodes = [m for m in self.massepunkte if m.active]
        current_count = len(active_nodes)
        target_count = int(len(self.massepunkte) * target_mass_percent)
        
        # Wie viele Massenpunkte müssen weg?
        to_remove_count = current_count - target_count
        
        if to_remove_count <= 0:
            return # Ziel schon erreicht
            
        candidates = []
        for m in active_nodes:
            # Knoten mit Kraft oder Randbedingung dürfen nicht gelöscht werden!
            is_fixed = np.any(m.fixed)
            has_force = np.linalg.norm(m.force) > 1e-9
            
            if not is_fixed and not has_force:
                candidates.append(m)
        
        # Aufsteigend sortieren nach Energie
        candidates.sort(key=lambda m: m.strain_energy)
        
        removed_count = 0
        
        for m in candidates:
            # Abbrechen, wenn wir genug gelöscht haben
            if removed_count >= to_remove_count:
                break 
                
            # 1. Testweise deaktivieren
            m.active = False
            
            # 2. Prüfen: Fällt das Fachwerk auseinander?
            if self.is_path_intact():
                # Pfad ist intakt! Löschung ist gültig.
                removed_count += 1
            else:
                # FEHLER! Pfad abgerissen. Sofort wiederherstellen!
                m.active = True
                m.strain_energy = 1e10 # Extrem hohe Energie geben, damit er nicht nochmal geprüft wird
                
        print(f"Optimierung: {removed_count} Knoten entfernt (Ziel für diesen Schritt war {to_remove_count}).")

        # Matrix wird neu aufgebaut, wenn das nächste Mal gelöst wird
        if hasattr(self, "K_base"):
            del self.K_base

    def generate_rect_mesh(self, nx: int, nz: int, width: float, height: float):
        """
        Erstellt ein Rechteck-Gitter mit nx * nz Knoten.
        """
        # 1. Alte Daten löschen
        self.massepunkte = []
        self.federn = []
        
        # Abstände berechnen
        dx = width / (nx - 1) if nx > 1 else 0
        dz = height / (nz - 1) if nz > 1 else 0
        
        # 2. Knoten erstellen
        for z_i in range(nz):
            for x_i in range(nx):
                # Koordinaten: x geht nach rechts, z nach unten
                self.add_Massepunkt(x_i * dx, z = z_i * dz)

        # 3. Federn erstellen
        # Wir gehen durch jeden Punkt und verbinden ihn mit seinen Nachbarn
        # (Rechts, Unten, Unten-Rechts, Unten-Links)
        
        k_ortho = 1.0
        k_diag = 1.0 / np.sqrt(2.0)
        
        for z_i in range(nz):
            for x_i in range(nx):
                # Aktueller Index im 1D-Array
                idx = z_i * nx + x_i
                
                # --- Nachbar Rechts (Horizontal) ---
                if x_i < nx - 1:
                    idx_right = idx + 1
                    self.add_Feder(idx, idx_right, k=k_ortho)
                
                # --- Nachbar Unten (Vertikal) ---
                if z_i < nz - 1:
                    idx_down = (z_i + 1) * nx + x_i
                    self.add_Feder(idx, idx_down, k=k_ortho)
                    
                # --- Nachbar Unten-Rechts (Diagonal) ---
                if x_i < nx - 1 and z_i < nz - 1:
                    idx_br = (z_i + 1) * nx + (x_i + 1)
                    self.add_Feder(idx, idx_br, k=k_diag)
                    
                # --- Nachbar Unten-Links (Diagonal) ---
                if x_i > 0 and z_i < nz - 1:
                    idx_bl = (z_i + 1) * nx + (x_i - 1)
                    self.add_Feder(idx, idx_bl, k=k_diag)
    
        # Matrix wird neu aufgebaut, wenn das nächste Mal Knoten entfernt oder gelöst wird
        if hasattr(self, "K_base"):
            del self.K_base
            
    def stable_test(self):

        #Freiheitsgrade und Lager überprüfen, damit die Optimierung durchgeführt werden kann.
        fixed_nodes = 0
        fixed_dofs = 0
        festlager = 0
        loslager = 0

        for m in self.massepunkte:
            if np.any(m.fixed):
                fixed_nodes += 1
                n_fixed = np.sum(m.fixed)
                fixed_dofs += n_fixed

                if n_fixed == self.dim:
                    festlager += 1
                else:
                    loslager += 1

        # Verschiedene Fehlermeldungen ausgeben, je nach Fehler
        # In 2D brauchen wir 3 fixierte Freiheitsgrade in 3D brauchen wir 6 fixierte Freiheitsgrade
        required_dofs = 3 if self.dim == 2 else 6

        if fixed_nodes < 2:
            return False, "Mindestens 2 Lager werden gebraucht!"

        if fixed_dofs < required_dofs:
            return False, f"Es werden mindestens {required_dofs} fixierte Freiheitsgrade benötigt!"
        if festlager < 1:
            return False, "Mindestens 1 Festlager wird benötigt!"
        
        
        return True, "Struktur stabil."


    def is_path_intact(self):
        """
        Prüft mit Graphentheorie (BFS), ob alle Lasten noch mit den Lagern verbunden sind.
        """
        import numpy as np
        
        active_nodes = [m for m in self.massepunkte if m.active]
        supports = [m for m in active_nodes if np.any(m.fixed)]
        forces = [m for m in active_nodes if np.linalg.norm(m.force) > 0]

        if not forces or not supports:
            return True # Nichts zu prüfen

        # 1. Adjazenzliste aufbauen, Verbindungen der aktiven Knoten durch aktive Federn
        adj = {id(m): [] for m in active_nodes}
        for f in self.federn:
            if f.massepunkt_i.active and f.massepunkt_j.active:
                adj[id(f.massepunkt_i)].append(id(f.massepunkt_j))
                adj[id(f.massepunkt_j)].append(id(f.massepunkt_i))

        # 2. Breitensuche von jedem Lager aus, um alle erreichbaren Knoten zu markieren
        visited = set(id(s) for s in supports)
        queue = [id(s) for s in supports]

        while queue:
            curr_id = queue.pop(0)
            for neighbor_id in adj.get(curr_id, []):
                if neighbor_id not in visited:
                    visited.add(neighbor_id)
                    queue.append(neighbor_id)

        # 3. Ergebnisprüfung: Sind die IDs aller Kraft-Knoten erreicht worden?
        for f_node in forces:
            if id(f_node) not in visited:
                return False # Pfad unterbrochen.
        
        return True
    

    def generate_box_mesh(self, nx: int, ny: int, nz: int, width: float, depth: float, height: float):
        """
        Erstellt ein 3D-Quader-Gitter mit nx * ny * nz Knoten.
        """
        # 1. Alte Daten löschen und Dimension auf 3D setzen
        self.massepunkte = []
        self.federn = []
        self.dim = 3
        
        # Abstände berechnen
        dx = width / (nx - 1) if nx > 1 else 0
        dy = depth / (ny - 1) if ny > 1 else 0
        dz = height / (nz - 1) if nz > 1 else 0
        
        # 2. Knoten erstellen
        for z_i in range(nz):
            for y_i in range(ny):
                for x_i in range(nx):
                    self.add_Massepunkt(x_i * dx, z_i * dz, y=y_i * dy)

        # Hilfsfunktion, 1D-Index zu 3D-Punkt berechnen
        def get_idx(x, y, z):
            return z * (nx * ny) + y * nx + x

        # 3. Federn erstellen
        directions = [
            (1, 0, 0), (0, 1, 0), (0, 0, 1),               # Orthogonal (Gerade)
            (1, 1, 0), (1, -1, 0),                         # Fläche XY (Diagonalen)
            (1, 0, 1), (-1, 0, 1),                         # Fläche XZ (Diagonalen)
            (0, 1, 1), (0, -1, 1),                         # Fläche YZ (Diagonalen)
            (1, 1, 1), (1, -1, 1), (-1, 1, 1), (-1, -1, 1) # Quer durch den Raum (Raumdiagonalen)
        ]
        
        for z_i in range(nz):
            for y_i in range(ny):
                for x_i in range(nx):
                    idx = get_idx(x_i, y_i, z_i)
                    
                    for dx_dir, dy_dir, dz_dir in directions:
                        nx_i = x_i + dx_dir
                        ny_i = y_i + dy_dir
                        nz_i = z_i + dz_dir
                        
                        # Prüfen, ob der Zielknoten noch im Gitter liegt
                        if 0 <= nx_i < nx and 0 <= ny_i < ny and 0 <= nz_i < nz:
                            idx_neighbor = get_idx(nx_i, ny_i, nz_i)
                            
                            # Steifigkeit abhängig von der Länge machen (wie 1/sqrt(2) im 2D)
                            length_factor = np.sqrt(dx_dir**2 + dy_dir**2 + dz_dir**2)
                            k_val = 1.0 / length_factor
                            
                            self.add_Feder(idx, idx_neighbor, k=k_val)
    
        # Wenn Matrix noch im Speicher ist, löschen um neuen Aufbau zu erzwingen
        if hasattr(self, "K_base"):
            del self.K_base