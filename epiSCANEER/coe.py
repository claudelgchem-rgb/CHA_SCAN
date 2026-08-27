import os, shutil, subprocess, sys
from . import lib_path

COE_DIR = os.path.join(lib_path, "coe")
CONF_NAME = "Energetics.properties"
CONF_PATH = os.path.join(COE_DIR, CONF_NAME)
# Generated (machine specific) copy of the config, kept out of the source tree
# so that the checked-in Energetics.properties is never modified.
GENERATED_CONF_DIR = os.path.join(COE_DIR, "generated_conf")
GENERATED_CONF_PATH = os.path.join(GENERATED_CONF_DIR, CONF_NAME)


def get_java_executable():
    """Locate a java runtime: $JAVA_HOME, then a bundled JRE, then $PATH."""
    java_home = os.environ.get("JAVA_HOME")
    if java_home:
        candidate = os.path.join(java_home, "bin", "java")
        if sys.platform == "win32":
            candidate += ".exe"
        if os.path.isfile(candidate):
            return candidate

    bundled = os.path.join(lib_path, "jre", "linux", "bin", "java")
    if os.path.isfile(bundled):
        return bundled

    found = shutil.which("java.exe" if sys.platform == "win32" else "java")
    if found:
        return found

    raise RuntimeError(
        "Java runtime not found. Install Java (>= 1.8) or set JAVA_HOME; "
        "the covariation step (covariance.algorithms.*) requires it."
    )


def check_conf_file():
    """Write a copy of Energetics.properties with HOME_DIRECTORY set to this checkout.

    The Java code loads Energetics.properties from the classpath and resolves
    its data files (data/Maxhom_McLachlan.metric) relative to HOME_DIRECTORY,
    which upstream ships as a path from the original author's machine.
    """
    with open(CONF_PATH) as f:
        line_list = f.readlines()

    for i, line in enumerate(line_list):
        if line.strip().startswith("HOME_DIRECTORY"):
            line_list[i] = "HOME_DIRECTORY=%s\n" % COE_DIR
            break
    else:
        line_list.append("HOME_DIRECTORY=%s\n" % COE_DIR)

    os.makedirs(GENERATED_CONF_DIR, exist_ok=True)
    with open(GENERATED_CONF_PATH, "w") as fo:
        fo.writelines(line_list)


class ProcCoe:
    # input aln, out coe file path are required
    def __init__(self, aln_file, coe_file, error_file, algorithm="McBASC"):
        self.algorithm = algorithm
        self.input_aln_file = os.path.abspath(aln_file)
        self.output_coe_file = os.path.abspath(coe_file)
        self.result = []

    # Available Algorithm: McBASC, ELSC, SCA, MI, OMES
    def run(self):
        check_conf_file()
        jexePath = get_java_executable()
        jClassPath = os.pathsep.join([GENERATED_CONF_DIR, COE_DIR])
        jClassAlgorithm = "covariance.algorithms.%s" % self.algorithm
        cmd = [jexePath, "-cp", jClassPath, jClassAlgorithm,
               self.input_aln_file, self.output_coe_file]
        proc = subprocess.run(cmd, cwd=COE_DIR, capture_output=True, text=True)
        if proc.returncode != 0 or not os.path.isfile(self.output_coe_file):
            raise RuntimeError(
                "covariation step failed (%s)\ncommand: %s\nstdout: %s\nstderr: %s"
                % (jClassAlgorithm, " ".join(cmd), proc.stdout, proc.stderr)
            )

    # Parse coe-calculation result
    def parse(self):
        f = open(self.output_coe_file, 'r')
        #f.next()
        f.readline()
        for line in f.readlines():
            fields = line.split()
            self.result.append((int(fields[0]), int(fields[1]), float(fields[2])))
        f.close()
        return self.result

    def convertResPos(self, res_pos_dic):
        f = open(self.output_coe_file, 'w')
        print("residue_1\tresidue_2\tscore", file=f)
        for cell in self.result:
            print('%d\t%d\t%s' %\
             (res_pos_dic[cell[0]]+1, res_pos_dic[cell[1]]+1, str(cell[2])), file=f)
        f.close()
