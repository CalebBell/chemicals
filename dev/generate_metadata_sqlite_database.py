import os
import sqlite3
import json
from pathlib import Path
from chemicals.identifiers import (
    load_chemical_preferences,
    ChemicalMetadataDB,
    folder,
    PUBCHEM_LARGE_DB_NAME,
    PUBCHEM_SMALL_DB_NAME,
    PUBCHEM_CATION_DB_NAME,
    PUBCHEM_ANION_DB_NAME,
    PUBCHEM_IONORGANIC_DB_NAME,
    PUBCHEM_EXAMPLE_DB_NAME
)
preferred = load_chemical_preferences()[0]


def initialize_db(db_path):
    """Create a new SQLite database with the chemical metadata schema"""
    
    # Delete existing database if it exists
    db_path = Path(db_path)
    if db_path.exists():
        db_path.unlink()
    
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    
    # Create tables with appropriate indexes
    cur.executescript("""
        -- Main chemicals table using CAS as primary key
        CREATE TABLE chemicals (
            cas INTEGER PRIMARY KEY,
            preferred BOOLEAN NOT NULL DEFAULT FALSE,
            pubchemid INTEGER,
            formula TEXT,
            mw REAL,
            smiles TEXT,
            inchi TEXT,
            inchi_key TEXT,
            iupac_name TEXT,
            common_name TEXT,
            raw_synonyms TEXT
        );
        CREATE INDEX formula_idx ON chemicals(formula);
        CREATE INDEX smiles_idx ON chemicals(smiles);
        CREATE INDEX inchi_idx ON chemicals(inchi);
        CREATE INDEX pubchemid_idx ON chemicals(pubchemid);
        CREATE INDEX inchi_key_idx ON chemicals(inchi_key);
        CREATE INDEX formula_preferred_idx ON chemicals(formula, preferred DESC);

        -- Normalized synonym lookup table
        CREATE TABLE chemical_synonyms (
            cas INTEGER,
            synonym TEXT,
            FOREIGN KEY (cas) REFERENCES chemicals(cas) ON DELETE CASCADE,
            PRIMARY KEY (cas, synonym)
        );
        
        CREATE INDEX synonym_idx ON chemical_synonyms(synonym);
    """)
    
    conn.commit()
    return conn

chemicals_added = set()
def add_chemical(cur, chemical):
    """Add a single ChemicalMetadata object to the database"""
    if chemical.CASs in chemicals_added:
        raise ValueError(f"Duplicate {chemical.CASs}")
    # Insert main chemical record
    cur.execute("""
        INSERT INTO chemicals (
            cas, preferred, pubchemid, formula, mw, smiles, inchi, 
            inchi_key, iupac_name, common_name, raw_synonyms
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        chemical.CAS,
        chemical.CASs in preferred,
        chemical.pubchemid,
        chemical.formula,
        chemical.MW,
        chemical.smiles,
        chemical.InChI,
        chemical.InChI_key,
        chemical.iupac_name,
        chemical.common_name,
        '\t'.join([v for v in chemical.synonyms[2:] if v])
    ))
    
    synonyms = [(chemical.CAS, syn) for syn in chemical.synonyms]
    for syn in chemical.synonyms:
        if syn.lower() != syn:
            synonyms.append((chemical.CAS, syn.lower()))
    
    cur.executemany(
        "INSERT OR IGNORE INTO chemical_synonyms (cas, synonym) VALUES (?, ?)",
        synonyms
    )
    chemicals_added.add(chemical.CASs)

def force_insert_synonyms(cur, chemical):
    """Force insert synonyms from a chemical, overwriting existing ones if necessary"""
    synonyms = [(chemical.CAS, syn) for syn in chemical.synonyms]
    for syn in chemical.synonyms:
        if syn.lower() != syn:
            synonyms.append((chemical.CAS, syn.lower()))
    
    for cas, syn in synonyms:
        cur.execute("DELETE FROM chemical_synonyms WHERE synonym = ?", (syn,))

    # Insert new synonyms
    cur.executemany(
        "INSERT OR REPLACE INTO chemical_synonyms (cas, synonym) VALUES (?, ?)",
        synonyms
    )

def dump_to_db(chemical_db, db_path='chemicals.db'):
    """Dump a ChemicalMetadataDB instance to SQLite database"""
    # small_dbs = ChemicalMetadataDB(elements=True, main_db=os.path.join(folder, PUBCHEM_SMALL_DB_NAME), 
    #                                 user_dbs=[os.path.join(folder, n) for n in [PUBCHEM_CATION_DB_NAME, PUBCHEM_ANION_DB_NAME,
    #                                         PUBCHEM_IONORGANIC_DB_NAME, PUBCHEM_EXAMPLE_DB_NAME]])
    small_dbs = ChemicalMetadataDB(elements=True, main_db=os.path.join(folder, PUBCHEM_LARGE_DB_NAME), 
                                    user_dbs=[os.path.join(folder, n) for n in [PUBCHEM_SMALL_DB_NAME, PUBCHEM_CATION_DB_NAME, PUBCHEM_ANION_DB_NAME,
                                            PUBCHEM_IONORGANIC_DB_NAME, PUBCHEM_EXAMPLE_DB_NAME]])
    conn = initialize_db(db_path)
    cur = conn.cursor()
    
    # Use a transaction for faster bulk insert
    cur.execute("BEGIN TRANSACTION")
    
    # try:
    # for chemical in chemical_db:
    for chemical in small_dbs:
        add_chemical(cur, chemical)
    # except Exception as e:
    #     conn.rollback()
    #     raise e
    # finally:
        # conn.close()

    # 2. Force insert synonyms from combined small DBs
    for chemical in small_dbs:
        force_insert_synonyms(cur, chemical)
    
    # 3. Force insert synonyms from other DBs only
    other_dbs = ChemicalMetadataDB(elements=True, main_db=None, 
                                    user_dbs=[os.path.join(folder, n) for n in [PUBCHEM_EXAMPLE_DB_NAME, PUBCHEM_IONORGANIC_DB_NAME, PUBCHEM_CATION_DB_NAME, PUBCHEM_ANION_DB_NAME]])
    for chemical in other_dbs:
        force_insert_synonyms(cur, chemical)
    conn.commit()
    conn.close()
        


def verify_db(db_path):
    """Print database statistics and verify integrity"""
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    
    print("Database Statistics:")
    
    # Check chemical counts
    cur.execute("SELECT COUNT(*) FROM chemicals")
    chemical_count = cur.fetchone()[0]
    print(f"Total chemicals: {chemical_count}")
    
    # Check synonym counts
    cur.execute("SELECT COUNT(*) FROM chemical_synonyms")
    synonym_count = cur.fetchone()[0]
    print(f"Total synonym entries: {synonym_count}")
    
    # Check for potential issues
    print("\nVerifying database integrity:")
    
    # Check for orphaned synonyms
    cur.execute("""
        SELECT COUNT(*) FROM chemical_synonyms cs 
        LEFT JOIN chemicals c ON cs.cas = c.cas 
        WHERE c.cas IS NULL
    """)
    orphans = cur.fetchone()[0]
    print(f"Orphaned synonyms: {orphans}")
    
    # Verify some key constraints
    cur.execute("PRAGMA integrity_check")
    integrity = cur.fetchone()[0]
    print(f"Database integrity: {integrity}")
    
    conn.close()

# Example usage:
if __name__ == '__main__':
    from chemicals.identifiers import pubchem_db, folder

    db_path = os.path.join(folder, 'metadata.db')
    pubchem_db.autoload_main_db()
    # Dump to SQLite
    dump_to_db(pubchem_db, db_path)
    
    # Verify the database
    verify_db(db_path)
    
    # Example of how to use the database:
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    
    # Search by name example
    name = "methanol"
    cur.execute("""
        SELECT c.* FROM chemicals c
        JOIN chemical_synonyms cs ON c.cas = cs.cas
        WHERE cs.synonym = ?
    """, (name,))
    result = cur.fetchone()
    
    # Search by CAS example
    cas = 67561  # Methanol's CAS without hyphens
    cur.execute("SELECT * FROM chemicals WHERE cas = ?", (cas,))
    result = cur.fetchone()
    print(result)
    
    conn.close()