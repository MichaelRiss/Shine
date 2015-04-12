/*  Shine - NMR NOESY simulation program
    Copyright (C) 2015  Michael Riss

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

#include <string.h>
#include <iostream>
#include "../include/amino.h"

using namespace std;

unordered_map<string, aminoflags> aflags;

void InitAminoFlags(){
  aflags["ALA"] = 
    ACYCLIC | ALIPHATIC | BURIED | HYDROPHOBIC | NEUTRAL | SMALL;
  aflags["ARG"] = 
    ACYCLIC | BASIC | CHARGED | LARGE | POLAR | POSITIVE | SURFACE;
  aflags["ASN"] = 
    ACYCLIC | MEDIUM | NEUTRAL | POLAR | SURFACE;
  aflags["ASP"] = 
    ACIDIC | ACYCLIC | CHARGED | MEDIUM | NEGATIVE | POLAR | SURFACE;
  aflags["CYS"] =
    ACYCLIC | BURIED | MEDIUM | NEUTRAL | POLAR;
  aflags["GLU"] = 
    ACIDIC | ACYCLIC | CHARGED | LARGE | NEGATIVE | POLAR | SURFACE;
  aflags["GLN"] =
    ACYCLIC | LARGE | NEUTRAL | POLAR | SURFACE;
  aflags["GLY"] =
    ACYCLIC | ALIPHATIC | HYDROPHOBIC | NEUTRAL | SMALL | SURFACE;
  aflags["HIS"] =
    AROMATIC | BASIC | CHARGED | CYCLIC | LARGE | NEUTRAL | POLAR | 
    POSITIVE | SURFACE;
  aflags["ILE"] =
    ACYCLIC | ALIPHATIC | BURIED | HYDROPHOBIC | LARGE | NEUTRAL;
  aflags["LEU"] =
    ACYCLIC | ALIPHATIC | BURIED | HYDROPHOBIC | LARGE | NEUTRAL;
  aflags["LYS"] =
    ACYCLIC | BASIC | CHARGED | LARGE | POLAR | POSITIVE | SURFACE;
  aflags["MET"] =
    ACYCLIC | BURIED | HYDROPHOBIC | LARGE | NEUTRAL;
  aflags["PHE"] =
    AROMATIC | BURIED | CYCLIC | HYDROPHOBIC | LARGE | NEUTRAL;
  aflags["PRO"] =
    CYCLIC | HYDROPHOBIC | MEDIUM | NEUTRAL | SURFACE;
  aflags["SER"] =
    ACYCLIC | NEUTRAL | POLAR | SMALL | SURFACE;
  aflags["THR"] =
    ACYCLIC | MEDIUM | NEUTRAL | POLAR | SURFACE;
  aflags["TRP"] =
    AROMATIC | BURIED | CYCLIC | HYDROPHOBIC | LARGE | NEUTRAL;
  aflags["TYR"] =
    AROMATIC | CYCLIC | HYDROPHOBIC | LARGE | NEUTRAL | SURFACE;
  aflags["VAL"] =
    ACYCLIC | ALIPHATIC | BURIED | HYDROPHOBIC | MEDIUM | NEUTRAL;
}

unordered_map<string, float> centroidsize;

void InitCentroidSize(){
  centroidsize["GLY"] = 1.0;
  centroidsize["ALA"] = 2.0;
  centroidsize["SER"] = 3.0;
  centroidsize["THR"] = 4.0;
  centroidsize["CYS"] = 5.0;
  centroidsize["VAL"] = 6.0;
  centroidsize["ASP"] = 7.0;
  centroidsize["PRO"] = 8.0;
  centroidsize["ILE"] = 9.0;
  centroidsize["LEU"] = 10.0;
  centroidsize["ASN"] = 11.0;
  centroidsize["PHE"] = 12.0;
  centroidsize["MET"] = 13.0;
  centroidsize["TYR"] = 14.0;
  centroidsize["GLU"] = 15.0;
  centroidsize["HIS"] = 16.0;
  centroidsize["TRP"] = 17.0;
  centroidsize["GLN"] = 18.0;
  centroidsize["LYS"] = 19.0;
  centroidsize["ARG"] = 20.0;
}

const char* contable[][9] = {
  { "A", "ALA", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "A", "ALA", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "A", "ALA", "HB1",  "1HB",  "MB",   "HB",   "CB",  "HB",   "QB"   },
  { "A", "ALA", "HB2",  "2HB",  "MB",   "HB",   "CB",  "HB",   "QB"   },
  { "A", "ALA", "HB3",  "3HB",  "MB",   "HB",   "CB",  "HB",   "QB"   },
  { "C", "CYS", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "C", "CYS", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "C", "CYS", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "C", "CYS", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "C", "CYS", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "C", "CYS", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "D", "ASP", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "D", "ASP", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "D", "ASP", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "D", "ASP", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "D", "ASP", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "D", "ASP", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "D", "ASP", "HB1",  "1HB",  "HBs",  "HB3",  "CB",  "HB",   "QB"   }, //
  { "D", "ASP", "HB2",  "2HB",  "HBs",  "HB2",  "CB",  "HB",   "QB"   }, //
  { "E", "GLU", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "E", "GLU", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "E", "GLU", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "E", "GLU", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "E", "GLU", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "E", "GLU", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "E", "GLU", "HB1",  "1HB",  "HBs",  "HB3",  "CB",  "HB",   "QB"   }, //
  { "E", "GLU", "HB2",  "2HB",  "HBs",  "HB2",  "CB",  "HB",   "QB"   }, //
  { "E", "GLU", "HG1",  "1HG",  "HGu",  "HG3",  "CG",  "HG",   "QG"   },
  { "E", "GLU", "HG2",  "2HG",  "HGd",  "HG2",  "CG",  "HG",   "QG"   },
  { "E", "GLU", "HG1",  "1HG",  "HG1",  "HG3",  "CG",  "HG1",  "HG3"  },
  { "E", "GLU", "HG2",  "2HG",  "HG2",  "HG2",  "CG",  "HG2",  "HG2"  },
  { "E", "GLU", "HG1",  "1HG",  "HGs",  "HG3",  "CG",  "HG",   "QG"   }, //
  { "E", "GLU", "HG2",  "2HG",  "HGs",  "HG2",  "CG",  "HG",   "QG"   }, //
  { "F", "PHE", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "F", "PHE", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "F", "PHE", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "F", "PHE", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "F", "PHE", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "F", "PHE", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "F", "PHE", "HD1",  "HD1",  "AD",   "HD1",  "CD1", "HD",   "QD"   },
  { "F", "PHE", "HD2",  "HD2",  "AD",   "HD2",  "CD2", "HD",   "QD"   },
  { "F", "PHE", "HE1",  "HE1",  "AE",   "HE1",  "CE1", "HE",   "QE"   },
  { "F", "PHE", "HE2",  "HE2",  "AE",   "HE2",  "CE2", "HE",   "QE"   },
  { "F", "PHE", "HZ",   "HZ",   "HZ",   "HZ",   "CZ",  "HZ",   "HZ"   },
  { "G", "GLY", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "G", "GLY", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "G", "GLY", "HA1",  "1HA",  "HAu",  "HA3",  "CA",  "HA",   "QA"   },
  { "G", "GLY", "HA2",  "2HA",  "HAd",  "HA2",  "CA",  "HA",   "QA"   },
  { "G", "GLY", "HA1",  "1HA",  "HA1",  "HA3",  "CA",  "HA1",  "HA1"  },
  { "G", "GLY", "HA2",  "2HA",  "HA2",  "HA2",  "CA",  "HA2",  "HA2"  },
  { "G", "GLY", "HA1",  "1HA",  "HAs",  "HA3",  "CA",  "HA",   "QA"   }, //
  { "G", "GLY", "HA2",  "2HA",  "HAs",  "HA2",  "CA",  "HA",   "QA"   }, //
  { "H", "HIS", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "H", "HIS", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "H", "HIS", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "H", "HIS", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "H", "HIS", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "H", "HIS", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "H", "HIS", "HD2",  "HD2",  "HD2",  "HD2",  "CD2", "HD2",  "HD2"  },
  { "H", "HIS", "HD1",  "HD1",  "HD1",  "HD1",  "ND1", "HD1",  "HD1"  },
  { "H", "HIS", "HE2",  "HE2",  "HE2",  "HE2",  "NE2", "HE2",  "HE2"  },
  { "H", "HIS", "HE1",  "HE1",  "HE1",  "HE1",  "CE1", "HE1",  "HE1"  },
  { "I", "ILE", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "I", "ILE", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "I", "ILE", "HB",   "HB",   "HB",   "HB",   "CB",  "HB",   "HB"   },
  { "I", "ILE", "HG11", "1HG1", "HGu",  "HG13", "CG1", "HG1",  "QG1"  },
  { "I", "ILE", "HG12", "2HG1", "HGd",  "HG12", "CG1", "HG1",  "QG1"  },
  { "I", "ILE", "HG11", "1HG1", "HG1",  "HG13", "CG1", "HG11", "HG13" },
  { "I", "ILE", "HG12", "2HG1", "HG2",  "HG12", "CG1", "HG12", "HG12" },
  { "I", "ILE", "HG11", "1HG1", "HGs",  "HG13", "CG1", "HG1",  "QG1"  }, //
  { "I", "ILE", "HG12", "2HG1", "HGs",  "HG12", "CG1", "HG1",  "QG1"  }, //
  { "I", "ILE", "HG21", "1HG2", "MG",   "HG2",  "CG2", "HG2",  "QG2"  },
  { "I", "ILE", "HG22", "2HG2", "MG",   "HG2",  "CG2", "HG2",  "QG2"  },
  { "I", "ILE", "HG23", "3HG2", "MG",   "HG2",  "CG2", "HG2",  "QG2"  },
  { "I", "ILE", "HD11", "1HD1", "MD",   "HD1",  "CD1", "HD1",  "QD1"  },
  { "I", "ILE", "HD12", "2HD1", "MD",   "HD1",  "CD1", "HD1",  "QD1"  },
  { "I", "ILE", "HD13", "3HD1", "MD",   "HD1",  "CD1", "HD1",  "QD1"  },
  { "K", "LYS", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "K", "LYS", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "K", "LYS", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "K", "LYS", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "K", "LYS", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "K", "LYS", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "K", "LYS", "HB1",  "1HB",  "HBs",  "HB3",  "CB",  "HB",   "QB"   }, //
  { "K", "LYS", "HB2",  "2HB",  "HBs",  "HB2",  "CB",  "HB",   "QB"   }, //
  { "K", "LYS", "HG1",  "1HG",  "HGu",  "HG3",  "CG",  "HG",   "QG"   },
  { "K", "LYS", "HG2",  "2HG",  "HGd",  "HG2",  "CG",  "HG",   "QG"   },
  { "K", "LYS", "HG1",  "1HG",  "HG1",  "HG3",  "CG",  "HG1",  "HG3"  },
  { "K", "LYS", "HG2",  "2HG",  "HG2",  "HG2",  "CG",  "HG2",  "HG2"  },
  { "K", "LYS", "HG1",  "1HG",  "HGs",  "HG3",  "CG",  "HG",   "QG"   }, //
  { "K", "LYS", "HG2",  "2HG",  "HGs",  "HG2",  "CG",  "HG",   "QG"   }, //
  { "K", "LYS", "HD1",  "1HD",  "HDu",  "HD3",  "CD",  "HD",   "QD"   },
  { "K", "LYS", "HD2",  "2HD",  "HDd",  "HD2",  "CD",  "HD",   "QD"   },
  { "K", "LYS", "HD1",  "1HD",  "HD1",  "HD3",  "CD",  "HD1",  "HD3"  },
  { "K", "LYS", "HD2",  "2HD",  "HD2",  "HD2",  "CD",  "HD2",  "HD2"  },
  { "K", "LYS", "HD1",  "1HD",  "HDs",  "HD3",  "CD",  "HD",   "QD"   }, //
  { "K", "LYS", "HD2",  "2HD",  "HDs",  "HD2",  "CD",  "HD",   "QD"   }, //
  { "K", "LYS", "HE1",  "1HE",  "HEu",  "HE3",  "CE",  "HE",   "QE"   },
  { "K", "LYS", "HE2",  "2HE",  "HEd",  "HE2",  "CE",  "HE",   "QE"   },
  { "K", "LYS", "HE1",  "1HE",  "HE1",  "HE3",  "CE",  "HE1",  "HE3"  },
  { "K", "LYS", "HE2",  "2HE",  "HE2",  "HE2",  "CE",  "HE2",  "HE2"  },
  { "K", "LYS", "HE1",  "1HE",  "HEs",  "HE3",  "CE",  "HE",   "QE"   }, //
  { "K", "LYS", "HE2",  "2HE",  "HEs",  "HE2",  "CE",  "HE",   "QE"   }, //
  { "K", "LYS", "HZ1",  "1HZ",  "HZ",   "HZ1",  "NZ",  "HZ",   "QZ"   },
  { "K", "LYS", "HZ2",  "2HZ",  "HZ",   "HZ2",  "NZ",  "HZ",   "QZ"   },
  { "K", "LYS", "HZ3",  "3HZ",  "HZ",   "HZ3",  "NZ",  "HZ",   "QZ"   },
  { "L", "LEU", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "L", "LEU", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "L", "LEU", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "L", "LEU", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "L", "LEU", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "L", "LEU", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "L", "LEU", "HG",   "HG",   "HG",   "HG",   "CG",  "HG",   "HG"   },
  { "L", "LEU", "HD11", "1HD1", "MDu",  "HD1",  "CDu", "HD",   "QD"   },
  { "L", "LEU", "HD12", "2HD1", "MDu",  "HD1",  "CDu", "HD",   "QD"   },
  { "L", "LEU", "HD13", "3HD1", "MDu",  "HD1",  "CDu", "HD",   "QD"   },
  { "L", "LEU", "HD21", "1HD2", "MDd",  "HD2",  "CDd", "HD",   "QD"   },
  { "L", "LEU", "HD22", "2HD2", "MDd",  "HD2",  "CDd", "HD",   "QD"   },
  { "L", "LEU", "HD23", "3HD2", "MDd",  "HD2",  "CDd", "HD",   "QD"   },
  { "L", "LEU", "HD11", "1HD1", "MD1",  "HD1",  "CD1", "HD1",  "QD1"  },
  { "L", "LEU", "HD12", "2HD1", "MD1",  "HD1",  "CD1", "HD1",  "QD1"  },
  { "L", "LEU", "HD13", "3HD1", "MD1",  "HD1",  "CD1", "HD1",  "QD1"  },
  { "L", "LEU", "HD21", "1HD2", "MD2",  "HD2",  "CD2", "HD2",  "QD2"  },
  { "L", "LEU", "HD22", "2HD2", "MD2",  "HD2",  "CD2", "HD2",  "QD2"  },
  { "L", "LEU", "HD23", "3HD2", "MD2",  "HD2",  "CD2", "HD2",  "QD2"  },
  { "L", "LEU", "HD11", "1HD1", "MDs",  "HD1",  "CDu", "HD",   "QD"   }, //
  { "L", "LEU", "HD12", "2HD1", "MDs",  "HD1",  "CDu", "HD",   "QD"   }, //
  { "L", "LEU", "HD13", "3HD1", "MDs",  "HD1",  "CDu", "HD",   "QD"   }, //
  { "L", "LEU", "HD21", "1HD2", "MDs",  "HD2",  "CDd", "HD",   "QD"   }, //
  { "L", "LEU", "HD22", "2HD2", "MDs",  "HD2",  "CDd", "HD",   "QD"   }, //
  { "L", "LEU", "HD23", "3HD2", "MDs",  "HD2",  "CDd", "HD",   "QD"   }, //
  { "M", "MET", "HT1",  "1HT",  "HT",   "H",    "NT",  "HT",   "HT"   },
  { "M", "MET", "HT2",  "2HT",  "HT",   "H",    "NT",  "HT",   "HT"   },
  { "M", "MET", "HT3",  "3HT",  "HT",   "H",    "NT",  "HT",   "HT"   },
  { "M", "MET", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "M", "MET", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "M", "MET", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "M", "MET", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "M", "MET", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "M", "MET", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "M", "MET", "HG1",  "1HG",  "HGu",  "HG3",  "CG",  "HG",   "QG"   },
  { "M", "MET", "HG2",  "2HG",  "HGd",  "HG2",  "CG",  "HG",   "QG"   },
  { "M", "MET", "HG1",  "1HG",  "HG1",  "HG3",  "CG",  "HG1",  "HG3"  },
  { "M", "MET", "HG2",  "2HG",  "HG2",  "HG2",  "CG",  "HG2",  "HG2"  },
  { "M", "MET", "HE1",  "1HE",  "ME",   "HE",   "CE",  "HE",   "QE"   },
  { "M", "MET", "HE2",  "2HE",  "ME",   "HE",   "CE",  "HE",   "QE"   },
  { "M", "MET", "HE3",  "3HE",  "ME",   "HE",   "CE",  "HE",   "QE"   },
  { "N", "ASN", "HT1",  "1HT",  "HT",   "H",    "NT",  "HT",   "HT"   },
  { "N", "ASN", "HT2",  "2HT",  "HT",   "H",    "NT",  "HT",   "HT"   },
  { "N", "ASN", "HT3",  "3HT",  "HT",   "H",    "NT",  "HT",   "HT"   },
  { "N", "ASN", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "N", "ASN", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "N", "ASN", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "N", "ASN", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "N", "ASN", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "N", "ASN", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "N", "ASN", "HD21", "1HD2", "HDu",  "HD21", "ND2", "HD2",  "QD2"  },
  { "N", "ASN", "HD22", "2HD2", "HDd",  "HD22", "ND2", "HD2",  "QD2"  },
  { "N", "ASN", "HD21", "1HD2", "HD1",  "HD21", "ND2", "HD21", "HD21" },
  { "N", "ASN", "HD22", "2HD2", "HD2",  "HD22", "ND2", "HD22", "HD22" },
  { "P", "PRO", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "P", "PRO", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "P", "PRO", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "P", "PRO", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "P", "PRO", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "P", "PRO", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "P", "PRO", "HG1",  "1HG",  "HGu",  "HG3",  "CG",  "HG",   "QG"   },
  { "P", "PRO", "HG2",  "2HG",  "HGd",  "HG2",  "CG",  "HG",   "QG"   },
  { "P", "PRO", "HG1",  "1HG",  "HG1",  "HG3",  "CG",  "HG1",  "HG3"  },
  { "P", "PRO", "HG2",  "2HG",  "HG2",  "HG2",  "CG",  "HG2",  "HG2"  },
  { "P", "PRO", "HG1",  "1HG",  "HGs",  "HG3",  "CG",  "HG",   "QG"   }, //
  { "P", "PRO", "HG2",  "2HG",  "HGs",  "HG2",  "CG",  "HG",   "QG"   }, //
  { "P", "PRO", "HD1",  "1HD",  "HDu",  "HD3",  "CD",  "HD",   "QD"   },
  { "P", "PRO", "HD2",  "2HD",  "HDd",  "HD2",  "CD",  "HD",   "QD"   },
  { "P", "PRO", "HD1",  "1HD",  "HD1",  "HD3",  "CD",  "HD1",  "HD3"  },
  { "P", "PRO", "HD2",  "2HD",  "HD2",  "HD2",  "CD",  "HD2",  "HD2"  },
  { "P", "PRO", "HD1",  "1HD",  "HDs",  "HD3",  "CD",  "HD",   "QD"   }, //
  { "P", "PRO", "HD2",  "2HD",  "HDs",  "HD2",  "CD",  "HD",   "QD"   }, //
  { "Q", "GLN", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "Q", "GLN", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "Q", "GLN", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "Q", "GLN", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "Q", "GLN", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "Q", "GLN", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "Q", "GLN", "HG1",  "1HG",  "HGu",  "HG3",  "CG",  "HG",   "QG"   },
  { "Q", "GLN", "HG2",  "2HG",  "HGd",  "HG2",  "CG",  "HG",   "QG"   },
  { "Q", "GLN", "HG1",  "1HG",  "HG1",  "HG3",  "CG",  "HG1",  "HG3"  },
  { "Q", "GLN", "HG2",  "2HG",  "HG2",  "HG2",  "CG",  "HG2",  "HG2"  },
  { "Q", "GLN", "HE21", "1HE2", "HEu",  "HE21", "NE2", "HE2",  "QE2"  },
  { "Q", "GLN", "HE22", "2HE2", "HEd",  "HE22", "NE2", "HE2",  "QE2"  },
  { "Q", "GLN", "HE21", "1HE2", "HE1",  "HE21", "NE2", "HE21", "HE21" },
  { "Q", "GLN", "HE22", "2HE2", "HE2",  "HE22", "NE2", "HE22", "HE22" },
  { "R", "ARG", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "R", "ARG", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "R", "ARG", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "R", "ARG", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "R", "ARG", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "R", "ARG", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "R", "ARG", "HB1",  "1HB",  "HBs",  "HB3",  "CB",  "HB",   "QB"   }, //
  { "R", "ARG", "HB2",  "2HB",  "HBs",  "HB2",  "CB",  "HB",   "QB"   }, //
  { "R", "ARG", "HG1",  "1HG",  "HGu",  "HG3",  "CG",  "HG",   "QG"   },
  { "R", "ARG", "HG2",  "2HG",  "HGd",  "HG2",  "CG",  "HG",   "QG"   },
  { "R", "ARG", "HG1",  "1HG",  "HG1",  "HG3",  "CG",  "HG1",  "HG3"  },
  { "R", "ARG", "HG2",  "2HG",  "HG2",  "HG2",  "CG",  "HG2",  "HG2"  },
  { "R", "ARG", "HD1",  "1HD",  "HDu",  "HD3",  "CD",  "HD",   "QD"   },
  { "R", "ARG", "HD2",  "2HD",  "HDd",  "HD2",  "CD",  "HD",   "QD"   },
  { "R", "ARG", "HD1",  "1HD",  "HD1",  "HD3",  "CD",  "HD1",  "HD3"  },
  { "R", "ARG", "HD2",  "2HD",  "HD2",  "HD2",  "CD",  "HD2",  "HD2"  },
  { "R", "ARG", "HD1",  "1HD",  "HDs",  "HD3",  "CD",  "HD",   "QD"   }, //
  { "R", "ARG", "HD2",  "2HD",  "HDs",  "HD2",  "CD",  "HD",   "QD"   }, //
  { "R", "ARG", "HE",   "HE",   "HE",   "HE",   "NE",  "HE",   "HE"   },
  { "R", "ARG", "HH11", "1HH1", "HH11", "HH11", "NH1", "HH1",  "QH1"  },
  { "R", "ARG", "HH12", "2HH1", "HH12", "HH12", "NH1", "HH1",  "QH1"  },
  { "R", "ARG", "NH1",  "NH1",  NULL,   "NH1",  "NH1", "NH1",  "NH1"  },
  { "R", "ARG", "HH21", "1HH2", "HH21", "HH21", "NH2", "HH2",  "QH2"  },
  { "R", "ARG", "HH22", "2HH2", "HH22", "HH22", "NH2", "HH2",  "QH2"  },
  { "R", "ARG", "NH2",  "NH2",  NULL,   "NH2",  "NH2", "NH2",  "NH2"  },
  { "S", "SER", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "S", "SER", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "S", "SER", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "S", "SER", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "S", "SER", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "S", "SER", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "S", "SER", "HG",   "HG",   "HG",   "HG",   "OG",  "HG",   "HG"   },
  { "T", "THR", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "T", "THR", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "T", "THR", "HB",   "HB",   "HB",   "HB",   "CB",  "HB",   "HB"   },
  { "T", "THR", "HG21", "1HG2", "MG2",  "HG2",  "CG2", "HG2",  "QG2"  },
  { "T", "THR", "HG22", "2HG2", "MG2",  "HG2",  "CG2", "HG2",  "QG2"  },
  { "T", "THR", "HG23", "3HG2", "MG2",  "HG2",  "CG2", "HG2",  "QG2"  },
  { "T", "THR", "HG1",  "1HG",  "HG1",  "HG11", "OG1", "HG1",  "HG11" },
  { "V", "VAL", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "V", "VAL", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "V", "VAL", "HB",   "HB",   "HB",   "HB",   "CB",  "HB",   "HB"   },
  { "V", "VAL", "HG11", "1HG1", "MGu",  "HG1",  "CGu", "HG",   "QG"   },
  { "V", "VAL", "HG12", "2HG1", "MGu",  "HG1",  "CGu", "HG",   "QG"   },
  { "V", "VAL", "HG13", "3HG1", "MGu",  "HG1",  "CGu", "HG",   "QG"   },
  { "V", "VAL", "HG21", "1HG2", "MGd",  "HG2",  "CGd", "HG",   "QG"   },
  { "V", "VAL", "HG22", "2HG2", "MGd",  "HG2",  "CGd", "HG",   "QG"   },
  { "V", "VAL", "HG23", "3HG2", "MGd",  "HG2",  "CGd", "HG",   "QG"   },
  { "V", "VAL", "HG11", "1HG1", "MG1",  "HG1",  "CG1", "HG1",  "QG1"  },
  { "V", "VAL", "HG12", "2HG1", "MG1",  "HG1",  "CG1", "HG1",  "QG1"  },
  { "V", "VAL", "HG13", "3HG1", "MG1",  "HG1",  "CG1", "HG1",  "QG1"  },
  { "V", "VAL", "HG21", "1HG2", "MG2",  "HG2",  "CG2", "HG2",  "QG2"  },
  { "V", "VAL", "HG22", "2HG2", "MG2",  "HG2",  "CG2", "HG2",  "QG2"  },
  { "V", "VAL", "HG23", "3HG2", "MG2",  "HG2",  "CG2", "HG2",  "QG2"  },
  { "W", "TRP", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "W", "TRP", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "W", "TRP", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "W", "TRP", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "W", "TRP", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "W", "TRP", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "W", "TRP", "CH2",  "CH2",  NULL,   "CH2",  "CH2", "CH2",  "CH2"  },
  { "W", "TRP", "HD1",  "HD1",  "HD1",  "HD1",  "CD1", "HD1",  "HD1"  },
  { "W", "TRP", "HE1",  "HE1",  "HE1",  "HE1",  "NE1", "HE1",  "HE1"  },
  { "W", "TRP", "HE3",  "HE3",  "HE3",  "HE3",  "CE3", "HE3",  "HE3"  },
  { "W", "TRP", "HZ2",  "HZ2",  "HZ2",  "HZ2",  "CZ2", "HZ2",  "HZ2"  },
  { "W", "TRP", "HZ3",  "HZ3",  "HZ3",  "HZ3",  "CZ3", "HZ3",  "HZ3"  },
  { "W", "TRP", "HH2",  "HH2",  "HH2",  "HH2",  "CH2", "HH2",  "HH2"  },
  { "Y", "TYR", "HN",   "H",    "HN",   "H",    "N",   "HN",   "HN"   },
  { "Y", "TYR", "HA",   "HA",   "HA",   "HA",   "CA",  "HA",   "HA"   },
  { "Y", "TYR", "HB1",  "1HB",  "HBu",  "HB3",  "CB",  "HB",   "QB"   },
  { "Y", "TYR", "HB2",  "2HB",  "HBd",  "HB2",  "CB",  "HB",   "QB"   },
  { "Y", "TYR", "HB1",  "1HB",  "HB1",  "HB3",  "CB",  "HB1",  "HB3"  },
  { "Y", "TYR", "HB2",  "2HB",  "HB2",  "HB2",  "CB",  "HB2",  "HB2"  },
  { "Y", "TYR", "HD1",  "HD1",  "AD",   "HD1",  "CD1", "HD",   "QD"   },
  { "Y", "TYR", "HD2",  "HD2",  "AD",   "HD2",  "CD2", "HD",   "QD"   },
  { "Y", "TYR", "HE1",  "HE1",  "AE",   "HE1",  "CE1", "HE",   "QE"   },
  { "Y", "TYR", "HE2",  "HE2",  "AE",   "HE2",  "CE2", "HE",   "QE"   },
  { "Y", "TYR", "OH",   "OH",   NULL,   "OH",   NULL,  "OH",   "OH"   },
  { "Y", "TYR", "HH",   "HH",   NULL,   NULL,   "OH",  NULL,   NULL   },
  { NULL, NULL, NULL,   NULL,   NULL,   NULL,   NULL,  NULL,   NULL   }
};

const char* methyltable[] = {
  "MB",
  "MG",
  "MD",
  "HZ",
  "MDu",
  "MDd",
  "MD1",
  "MD2",
  "MDs",
  "HT",
  "ME",
  "MG1",
  "MG2",
  "MGu",
  "MGd",
  NULL
};

const char* ac_oneletter_longlist2heavy( const char oneletter, 
					 const char* longlist ){
  int i = 0;
  while( contable[i][0] != NULL ){
    if( oneletter == *contable[i][0] && 
	strncmp( longlist, contable[i][4], 3 ) == 0 ){
      return contable[i][6];
    }
    i++;
  }
  return NULL;
}

const char* ac_oneletter_xplor2heavy( const char oneletter, 
				      const char* xplor ){
  int i = 0;
  while( contable[i][0] != NULL ){
    if( oneletter == *contable[i][0] && 
	strncmp( xplor, contable[i][2], 4 ) == 0 ){
      return contable[i][6];
    }
    i++;
  }
  return NULL;
}

const char* ac_oneletter2threeletter( const char oneletter ){
  int i = 0;
  while( contable[i][0] != NULL ){
    if( oneletter == *contable[i][0] ){
      return contable[i][1];
    }
    i++;
  }
  return NULL;
}

const char ac_threeletter2oneletter( const char* threeletter ){
  int i = 0;
  while( contable[i][0] != NULL ){
    if( strncmp( threeletter, contable[i][1], 3 ) == 0 ){
      return *contable[i][0];
    }
    i++;
  }
  return ' ';
}

const vector<const char*> 
ac_threeletter_xplor2longlist( const char* threeletter,
			       const char* xplor ){
  int i = 0;
  vector<const char*> result;
  while( contable[i][0] != NULL ){
    if( strncmp( threeletter, contable[i][1], 3 ) == 0 &&
	strncmp( xplor, contable[i][2], 4 ) == 0 ){
      result.push_back( contable[i][4] );
    }
    i++;
  }
  return result;
}

const char* ac_threeletter_xplor2heavy( const char* threeletter,
					const char* xplor ){
  int i = 0;
  while( contable[i][0] != NULL ){
    cout << "Compare " << threeletter << " with "
 	 << contable[i][1] << " and " << xplor
 	 << " with " << contable[i][2] << endl;
    if( strncmp( threeletter, contable[i][1], 3 ) == 0 &&
	strncmp( xplor, contable[i][2], 4 ) == 0 ){
      return contable[i][6];
    }
    i++;
  }
  cout << "Oops, returning NULL." << endl;
  return NULL;
}

const char* ac_threeletter_iupac2xplor( const char* threeletter,
					const char* iupac ){
  int i = 0;
  while( contable[i][0] != NULL ){
    if( strncmp( threeletter, contable[i][1], 3 ) == 0 &&
	strncmp( iupac, contable[i][3], 4 ) == 0 ){
      return contable[i][2];
    }
    i++;
  }
  cout << "ac_iupac2xplor: Oops, returning NULL." << endl;
  return NULL;
}

bool isMethylAtom( const char* atom ){
  // better use longlist format or disaster happens ...
  int i = 0;
  while( methyltable[i] != NULL ){
    if( strncmp( atom, methyltable[i], 3 ) ){
      return true;
    }
    i++;
  }
  return false;
}
