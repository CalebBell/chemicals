import sqlite3
import json
from pathlib import Path

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
            pubchemid INTEGER,
            formula TEXT,
            mw REAL,
            smiles TEXT,
            inchi TEXT,
            inchi_key TEXT,
            iupac_name TEXT,
            common_name TEXT,
            raw_synonyms TEXT  -- JSON array of original synonyms
        );
        CREATE INDEX formula_idx ON chemicals(formula);
        CREATE INDEX smiles_idx ON chemicals(smiles);
        CREATE INDEX inchi_idx ON chemicals(inchi);
        CREATE INDEX pubchemid_idx ON chemicals(pubchemid);
        CREATE INDEX inchi_key_idx ON chemicals(inchi_key);
        
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

def add_chemical(cur, chemical):
    """Add a single ChemicalMetadata object to the database"""
    
    # Insert main chemical record
    cur.execute("""
        INSERT INTO chemicals (
            cas, pubchemid, formula, mw, smiles, inchi, 
            inchi_key, iupac_name, common_name, raw_synonyms
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        chemical.CAS,
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

def dump_to_db(chemical_db, db_path='chemicals.db'):
    """Dump a ChemicalMetadataDB instance to SQLite database"""
    conn = initialize_db(db_path)
    cur = conn.cursor()
    
    # Use a transaction for faster bulk insert
    cur.execute("BEGIN TRANSACTION")
    
    try:
        for chemical in chemical_db:
            add_chemical(cur, chemical)
        conn.commit()
    except Exception as e:
        conn.rollback()
        raise e
    finally:
        conn.close()

def verify_db(db_path='chemicals.db'):
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
    from chemicals.identifiers import pubchem_db

    db_path = '../chemicals/Identifiers/metadata.db'
    pubchem_db.autoload_main_db()
    # Dump to SQLite
    dump_to_db(pubchem_db, db_path)
    
    # Verify the database
    verify_db()
    
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