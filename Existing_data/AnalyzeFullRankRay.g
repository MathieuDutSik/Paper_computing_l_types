# Analysis of the candidate rigid Gram matrices produced by
# LATT_LookForFullRankRayDomain and stored in DATA_LookForFullRankRay_dim6.
#
# For every matrix in the directory the script:
#   1. Computes its GL_n(Z) canonical form with LATT_Canonicalize and groups
#      the matrices into arithmetic-equivalence classes (one representative
#      per class).
#   2. Builds the iso-Delaunay domain of each representative with
#      LATT_SerialComputeDelaunay (QUERIES option FileIsoDelaunayDomain, which
#      dumps the domain as a boost archive).
#   3. Reads back that archive with LATT_AnalysisIsoDelaunay and records the
#      number of defining inequalities (full set n_ineq and irredundant/facet
#      set n_irred) together with the number of extreme rays n_ray.
#
# The three programs LATT_Canonicalize, LATT_SerialComputeDelaunay and
# LATT_AnalysisIsoDelaunay are taken from the PATH; it is the user's
# responsibility to make them available there (e.g. by adding the relevant
# polyhedral_common src_* directories to PATH).
#
# Run it from this directory:
#     gap AnalyzeFullRankRay.g
#
# POLY_MAX_CASES (default 0 = all) caps the number of representatives that go
# through the (expensive) steps 2-3, which is handy for a quick partial run.

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
data_dir := "DATA_LookForFullRankRay_dim6";

prog_canon := "LATT_Canonicalize";
prog_delaunay := "LATT_SerialComputeDelaunay";
prog_ana := "LATT_AnalysisIsoDelaunay";

max_cases := 0;
if IsBound(GAPInfo.SystemEnvironment.POLY_MAX_CASES) then
    max_cases := Int(GAPInfo.SystemEnvironment.POLY_MAX_CASES);
fi;

TmpDir := DirectoryTemporary();
tmp := function(name)
    return Filename(TmpDir, name);
end;

remove_if_exists := function(eFile)
    if IsExistingFile(eFile) then
        RemoveFile(eFile);
    fi;
end;

# The programs are looked up in the PATH via "which"; error out early with a
# clear message if one of them is missing.
CheckProgram := function(prog)
    local FileTmp;
    FileTmp := tmp("which.out");
    remove_if_exists(FileTmp);
    Exec(Concatenation("which ", prog, " > ", FileTmp, " 2> /dev/null"));
    if IsExistingFile(FileTmp) = false or Length(StringFile(FileTmp)) = 0 then
        Print("The program ", prog, " is not in the PATH.\n");
        Print("Make it available (e.g. add its polyhedral_common src_* ",
              "directory to PATH).\n");
        Error("Missing binary");
    fi;
    remove_if_exists(FileTmp);
end;
CheckProgram(prog_canon);
CheckProgram(prog_delaunay);
CheckProgram(prog_ana);

# --------------------------------------------------------------------------- #
# Step 1: arithmetic-equivalence classes via LATT_Canonicalize
# --------------------------------------------------------------------------- #

# GL_n(Z) canonical form of the Gram matrix stored in GramFile, returned as a
# GAP matrix. Two inputs are arithmetically equivalent iff these agree.
canonical_form := function(GramFile)
    local FileO, command;
    FileO := tmp("canon.out");
    remove_if_exists(FileO);
    command:=Concatenation(prog_canon, " gmp ", GramFile, " GAP ", FileO);
    Print("command=", command, "\n");
    Exec(command);
    if IsExistingFile(FileO) = false then
        Print("LATT_Canonicalize produced no output for ", GramFile, "\n");
        Error("Canonicalization failed");
    fi;
    return ReadAsFunction(FileO)();
end;

# The regular files of data_dir (everything except "." and ".."), sorted so
# that the class ordering and the chosen representatives are reproducible.
list_data_files := function()
    return SortedList(Filtered(DirectoryContents(data_dir),
                               x -> x <> "." and x <> ".."));
end;

# Group the files of data_dir by canonical form. Each class records the
# representative file, its canonical form and the list of member files.
compute_classes := function()
    local classes, eFile, GramFile, C, pos, i;
    classes := [];
    for eFile in list_data_files() do
        GramFile := Concatenation(data_dir, "/", eFile);
        C := canonical_form(GramFile);
        pos := fail;
        for i in [1..Length(classes)] do
            if classes[i].canon = C then
                pos := i;
                break;
            fi;
        od;
        if pos = fail then
            Add(classes, rec(canon := C, rep := eFile, members := [eFile]));
        else
            Add(classes[pos].members, eFile);
        fi;
    od;
    return classes;
end;

# --------------------------------------------------------------------------- #
# Step 2: build the iso-Delaunay domain with LATT_SerialComputeDelaunay
# --------------------------------------------------------------------------- #

# Runs LATT_SerialComputeDelaunay on GramFile and dumps the iso-Delaunay domain
# to ArchiveFile (QUERIES option FileIsoDelaunayDomain). Returns whether the
# archive was produced (it is not when the form is not fully rigid, i.e. when
# rigidity <> n*(n+1)/2).
build_iso_delaunay := function(GramFile, ArchiveFile)
    local FileNml, FileErr, command, os;
    FileNml := tmp("delaunay.nml");
    FileErr := tmp("delaunay.err");
    remove_if_exists(FileNml);
    remove_if_exists(ArchiveFile);
    os := OutputTextFile(FileNml, false);
    # Long temp-dir paths must not be wrapped in the middle of a "..." field.
    SetPrintFormattingStatus(os, false);
    AppendTo(os, "&SYSTEM\n");
    AppendTo(os, " OutFormat = \"nothing\"\n");
    AppendTo(os, " OutFile = \"", tmp("delaunay.unset.out"), "\"\n");
    AppendTo(os, " max_runtime_second = 0\n");
    AppendTo(os, "/\n\n");
    AppendTo(os, "&DATA\n");
    AppendTo(os, " arithmetic = \"gmp\"\n");
    AppendTo(os, " GRAMfile = \"", GramFile, "\"\n");
    AppendTo(os, " SVRfile = \"unset.svr\"\n");
    AppendTo(os, " CacheFile = \"none\"\n");
    AppendTo(os, "/\n\n");
    AppendTo(os, "&QUERIES\n");
    AppendTo(os, " FileIsoDelaunayDomain = \"", ArchiveFile, "\"\n");
    AppendTo(os, "/\n");
    CloseStream(os);
    command:=Concatenation(prog_delaunay, " ", FileNml, " 2> ", FileErr);
    Print("command=", command, "\n");
    Exec(command);
    remove_if_exists(FileNml);
    return IsExistingFile(ArchiveFile);
end;

# --------------------------------------------------------------------------- #
# Step 3: statistics with LATT_AnalysisIsoDelaunay
# --------------------------------------------------------------------------- #

# Parse the "key=value" report written by LATT_AnalysisIsoDelaunay into a record
# (fields n, sym_dim, n_ineq, n_irred, n_ray as integers).
parse_analysis := function(FileOut)
    local res, tok, parts;
    res := rec();
    for tok in SplitString(StringFile(FileOut), " \n") do
        parts := SplitString(tok, "=");
        if Length(parts) = 2 and Length(parts[1]) > 0 then
            res.(parts[1]) := Int(parts[2]);
        fi;
    od;
    return res;
end;

# Runs LATT_AnalysisIsoDelaunay on the iso-Delaunay archive. The FileRankTally
# query forces the (potentially expensive) extreme-ray enumeration, so that the
# report also contains n_ray. Returns the parsed record.
analyze_iso_delaunay := function(ArchiveFile)
    local FileNml, FileOut, FileTally, FileErr, command, os;
    FileNml := tmp("ana.nml");
    FileOut := tmp("ana.out");
    FileTally := tmp("ana.tally");
    FileErr := tmp("ana.err");
    remove_if_exists(FileNml);
    remove_if_exists(FileOut);
    os := OutputTextFile(FileNml, false);
    SetPrintFormattingStatus(os, false);
    AppendTo(os, "&DATA\n");
    AppendTo(os, " arithmetic = \"gmp\"\n");
    AppendTo(os, " FileIsoDelaunay = \"", ArchiveFile, "\"\n");
    AppendTo(os, "/\n\n");
    AppendTo(os, "&SYSTEM\n");
    AppendTo(os, " OutFile = \"", FileOut, "\"\n");
    AppendTo(os, "/\n\n");
    AppendTo(os, "&QUERIES\n");
    AppendTo(os, " FileRankTally = \"", FileTally, "\"\n");
    AppendTo(os, "/\n");
    CloseStream(os);
    command:=Concatenation(prog_ana, " ", FileNml, " 2> ", FileErr);
    Print("command=", command, "\n");
    Exec(command);
    remove_if_exists(FileNml);
    if IsExistingFile(FileOut) = false then
        Print("LATT_AnalysisIsoDelaunay produced no output for ", ArchiveFile,
              "\n");
        Error("Analysis failed");
    fi;
    return parse_analysis(FileOut);
end;

# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #

# Save a GAP object as "return OBJ;" to FileName.
save_data := function(FileName, OBJ)
    local os;
    remove_if_exists(FileName);
    os := OutputTextFile(FileName, false);
    AppendTo(os, "return ", OBJ, ";\n");
    CloseStream(os);
end;

FullRun := function()
    local classes, n_classes, n_todo, results, i, eClass, GramFile, Archive,
          ok, ana, eEntry;
    Print("Reading matrices from ", data_dir, "\n");
    classes := compute_classes();
    n_classes := Length(classes);
    Print("Number of matrices grouped into ", n_classes,
          " arithmetic-equivalence class(es)\n");
    Print("Class sizes: ", List(classes, x -> Length(x.members)), "\n\n");

    n_todo := n_classes;
    if max_cases > 0 and max_cases < n_todo then
        n_todo := max_cases;
        Print("POLY_MAX_CASES=", max_cases, ": only the first ", n_todo,
              " representative(s) are processed\n\n");
    fi;

    results := [];
    for i in [1..n_todo] do
        eClass := classes[i];
        GramFile := Concatenation(data_dir, "/", eClass.rep);
        Archive := tmp(Concatenation("iso_", String(i)));
        Print("i=", i, "/", n_todo, " rep=", eClass.rep,
              " (class size ", Length(eClass.members), ")\n");
        ok := build_iso_delaunay(GramFile, Archive);
        if ok = false then
            Print("  iso-Delaunay domain not built (form is not fully ",
                  "rigid?), skipping\n");
            eEntry := rec(rep := eClass.rep, class_size := Length(eClass.members),
                          members := eClass.members, built := false);
        else
            ana := analyze_iso_delaunay(Archive);
            remove_if_exists(Archive);
            Print("  n_ineq=", ana.n_ineq, " n_irred=", ana.n_irred,
                  " n_ray=", ana.n_ray, "\n");
            eEntry := rec(rep := eClass.rep, class_size := Length(eClass.members),
                          members := eClass.members, built := true,
                          n := ana.n, sym_dim := ana.sym_dim,
                          n_ineq := ana.n_ineq, n_irred := ana.n_irred,
                          n_ray := ana.n_ray);
        fi;
        Add(results, eEntry);
    od;
    return rec(n_classes := n_classes, classes := classes, results := results);
end;

Summary := function(full)
    local built, n_irred_list, n_ray_list;
    Print("\n=========================== SUMMARY ===========================\n");
    Print("Equivalence classes: ", full.n_classes, "\n");
    built := Filtered(full.results, x -> x.built);
    Print("Iso-Delaunay domains built and analyzed: ", Length(built), "\n");
    if Length(built) > 0 then
        n_irred_list := List(built, x -> x.n_irred);
        n_ray_list := List(built, x -> x.n_ray);
        Print("Defining inequalities (irredundant/facets) n_irred: ",
              "min=", Minimum(n_irred_list), " max=", Maximum(n_irred_list),
              " distinct=", Set(n_irred_list), "\n");
        Print("Extreme rays n_ray: min=", Minimum(n_ray_list),
              " max=", Maximum(n_ray_list), " distinct=", Set(n_ray_list),
              "\n");
    fi;
    Print("Detailed results saved to Result_FullRankRay.g\n");
    Print("===============================================================\n");
end;

full := FullRun();
Summary(full);
save_data("Result_FullRankRay.g", full.results);
