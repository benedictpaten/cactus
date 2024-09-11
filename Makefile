rootPath = .

include ${rootPath}/include.mk

modules = api setup caf bar hal reference pipeline preprocessor

# submodules are in multiple pass to handle dependencies cactus2hal being dependent on
# both cactus and sonLib
submodules1 = sonLib cPecan hal matchingAndOrdering pinchesAndCacti abPOA lastz paffy red FASTGA
submodules2 = cactus2hal
submodules = ${submodules1} ${submodules2}

git_commit ?= $(shell git rev-parse HEAD)
dockstore = quay.io/comparative-genomics-toolkit
name = ${dockstore}/cactus
tag = ${git_commit}
runtime_fullpath = $(realpath runtime)

# FIXME: this is also built by setup.py need to sort this out
versionPy = src/cactus/shared/version.py

CWD = ${PWD}

# these must be absolute, as used in submodules.
export sonLibRootDir = ${CWD}/submodules/sonLib
.PHONY: all all.% clean clean.% selfClean suball suball.% subclean.%

##
# Building.  First build submodules, then a pass for libs and a pass for bins
##
all:
	${MAKE} suball1
	${MAKE} all_libs
	${MAKE} all_progs
	${MAKE} suball2

# Note: hdf5 from apt doesn't seem to work for static builds.  It should be installed
# from source and configured with "--enable-static --disable-shared", then have its
# bin put at the front of PATH
static:
	CFLAGS="$${CFLAGS} -static" \
	CXXFLAGS="$${CXXFLAGS} -static" \
	LIBS="-static" \
	KC_OPTIONS="--enable-lzo --enable-static --disable-shared" \
	KT_OPTIONS="--without-lua --enable-static --disable-shared" \
	${MAKE} all

check-static: static
	if [ $(shell ls bin/* | grep -v mash | xargs ldd 2>& 1 | grep "not a dynamic" | wc -l) = $(shell ls bin/* | grep -v mash | wc -l) ] ; then\
		echo "ldd verified that all files in bin/ are static";\
	else\
		echo "ldd found dynamic linked binary in bin/";\
		exit 1;\
	fi

all_libs:
	${MAKE} ${modules:%=all_libs.%}
all_progs: all_libs
	${MAKE} ${modules:%=all_progs.%}

all_libs.%:
	cd $* && ${MAKE} all_libs
all_progs.%:
	cd $* && ${MAKE} all_progs

# special ordering
all_libs.blastLib: all_libs.api

##
# tests, see DEVELOPMENT.md for environment variables controling tests.
##

# Python tests
testModules = \
    progressive/outgroupTest.py \
    preprocessor/cactus_preprocessorTest.py \
    preprocessor/lastzRepeatMasking/cactus_lastzRepeatMaskTest.py \
    progressive/multiCactusTreeTest.py

# Unit tests (just collecting everything in bin/ with "test" in the name)
unitTests = \
	cactus_barTests \
	cactusAPITests \
	cactus_halGeneratorTests \
	stCafTests \
	stPinchesAndCactiTests \
	stPipelineTests \
	matchingAndOrderingTests \
	referenceTests \
	cPecanLibTests \
	sonLibTests \

#paffyTests \ # This is removed for now

# if running travis or gitlab, we want output to go to stdout/stderr so it can
# be seen in the log file, as opposed to individual files, which are much
# easier to read when running tests in parallel.
ifeq (${TRAVIS}${GITLAB_CI},) 
define testErrOut
     >& ${testLogDir}/$(basename $(notdir $*)).log
endef
else
    testErrOut =
endif

export PATH := ${CWD}/bin:${PATH}
export PYTHONPATH = ${CWD}/submodules/sonLib/src:${CWD}/submodules:${CWD}/src
pytestOpts = --tb=native --durations=0 -rsx
testOutDir = test-output
testLogDir = ${testOutDir}/logs

# parallel tests don't currenty work, not all cases of collission have
# been fixed
.NOTPARALLEL: test test_blast test_nonblast
test: ${testModules:%=%_runtest} ${unitTests:%=%_run_unit_test}
test_blast: ${testModules:%=%_runtest_blast}
test_nonblast: ${testModules:%=%_runtest_nonblast}
hal_test:
	cd ${CWD}/submodules/hal && make test

# run one test and save output
%_runtest: ${versionPy}
	@mkdir -p ${testLogDir}
	${PYTHON} -m pytest ${pytestOpts} src/cactus/$* ${testErrOut}

%_runtest_blast: ${versionPy}
	@mkdir -p ${testLogDir}
	${PYTHON} -m pytest ${pytestOpts} src/cactus/$* --suite=blast ${testErrOut} 

%_runtest_nonblast: ${versionPy}
	@mkdir -p ${testLogDir}
	${PYTHON} -m pytest ${pytestOpts} src/cactus/$* --suite=nonblast ${testErrOut}

%_run_unit_test:
	$*

${versionPy}:
	echo "cactus_commit = '${git_commit}'" >$@

${CWD}/test/mammals-truth.maf:
	cd ${CWD}/test && wget -q https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/all.maf -O mammals-truth.maf

${CWD}/test/primates-truth.maf:
	cd ${CWD}/test && wget -q https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/primates/loci1/all.maf -O primates-truth.maf

evolver_test: all ${CWD}/test/mammals-truth.maf ${CWD}/test/primates-truth.maf
# note: make docker needs to be run beforehand
	PYTHONPATH="${CWD}/submodules/" CACTUS_DOCKER_ORG=evolvertestdocker CACTUS_USE_LATEST=1 ${PYTHON} -m pytest ${pytestOpts} test/evolverTest.py

evolver_test_local: all ${CWD}/test/mammals-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverLocal

evolver_test_prepare_wdl: all ${CWD}/test/mammals-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverPrepareWDL

evolver_test_prepare_toil: all ${CWD}/test/mammals-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverPrepareToil

evolver_test_decomposed_local: all ${CWD}/test/mammals-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverDecomposedLocal

evolver_test_decomposed_docker: all ${CWD}/test/mammals-truth.maf
#note make docker needs to be run beforehand
	PYTHONPATH="${CWD}/submodules/" CACTUS_DOCKER_ORG=evolvertestdocker CACTUS_USE_LATEST=1 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverDecomposedDocker

evolver_test_docker: all ${CWD}/test/mammals-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverDocker

evolver_test_prepare_no_outgroup_docker: all ${CWD}/test/mammals-truth.maf
#note make docker needs to be run beforehand
	PYTHONPATH="${CWD}/submodules/" CACTUS_DOCKER_ORG=evolvertestdocker CACTUS_USE_LATEST=1 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverPrepareNoOutgroupDocker

evolver_test_prepare_no_outgroup_local: all ${CWD}/test/mammals-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverPrepareNoOutgroupLocal

evolver_test_update_node_local: ${CWD}/test/mammals-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverUpdateNodeLocal

evolver_test_update_branch_local: ${CWD}/test/mammals-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverUpdateBranchLocal

evolver_test_poa_local: all ${CWD}/test/primates-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverPOALocal

evolver_test_refmap_local: all ${CWD}/test/primates-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverRefmapLocal

evolver_test_minigraph_local: all ${CWD}/test/primates-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverMinigraphLocal

evolver_test_primates_pangenome_local: all ${CWD}/test/primates-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverPrimatesPangenomeLocal

evolver_test_primates_pangenome_docker: all ${CWD}/test/primates-truth.maf
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=docker CACTUS_DOCKER_MODE=1 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testEvolverPrimatesPangenomeDocker


evolver_test_all_local: evolver_test_local evolver_test_prepare_toil evolver_test_decomposed_local evolver_test_prepare_no_outgroup_local evolver_test_poa_local evolver_test_refmap_local evolver_test_minigraph_local

yeast_test_local:
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testYeastPangenomeLocal

yeast_test_step_by_step_local:
	PYTHONPATH="${CWD}/submodules/" CACTUS_BINARIES_MODE=local CACTUS_DOCKER_MODE=0 ${PYTHON} -m pytest ${pytestOpts} -s test/evolverTest.py::TestCase::testYeastPangenomeStepByStepLocal

##
# clean targets
##
selfClean: ${modules:%=clean.%}
	rm -rf include lib bin libexec ${versionPy} ${testOutDir} testTCDatabase build

clean.%:
	cd $* && ${MAKE} clean

clean: selfClean ${submodules:%=subclean.%}

##
# submodules
##
suball1: ${submodules1:%=suball.%}
suball2: ${submodules2:%=suball.%}

suball.sonLib:
	cd submodules/sonLib && PKG_CONFIG_PATH=${CWD}/lib/pkgconfig:${PKG_CONFIG_PATH} ${MAKE}
	mkdir -p ${BINDIR} ${LIBDIR} ${INCLDIR}
	rm -rf submodules/sonLib/bin/*.dSYM
	ln -f submodules/sonLib/bin/[a-zA-Z]* ${BINDIR}
	ln -f submodules/sonLib/lib/*.a ${LIBDIR}

suball.pinchesAndCacti: suball.sonLib
	cd submodules/pinchesAndCacti && ${MAKE}

suball.matchingAndOrdering: suball.sonLib
	cd submodules/matchingAndOrdering && ${MAKE}

suball.cPecan: suball.sonLib
	cd submodules/cPecan && ${MAKE}
	rm -f ${BINDIR}/cPecanLastz*

suball.cactus2hal: suball.sonLib suball.hal all_libs.api
	cd submodules/cactus2hal && ${MAKE}
	mkdir -p bin
	-ln -f submodules/cactus2hal/bin/* bin/

suball.hal: suball.sonLib
	cd submodules/hal &&  ${MAKE}
	mkdir -p bin
	-ln -f submodules/hal/bin/* bin/
	-ln -f submodules/hal/lib/libHal.a submodules/hal/lib/halLib.a

suball.abPOA:
	cd submodules/abPOA && ${MAKE}
	ln -f submodules/abPOA/lib/*.a ${LIBDIR}
	ln -f submodules/abPOA/include/*.h ${INCLDIR}
	rm -fr ${INCLDIR}/simde && cp -r submodules/abPOA/include/simde ${INCLDIR}

suball.lastz:
	cd submodules/lastz && ${MAKE}
	mkdir -p bin
	ln -f submodules/lastz/src/lastz bin

suball.paffy:
	cd submodules/paffy && ${MAKE}
	rm -rf submodules/paffy/bin/*.dSYM
	mkdir -p ${BINDIR} ${LIBDIR} ${INCLDIR}
	ln -f submodules/paffy/bin/[a-zA-Z]* ${BINDIR}
	ln -f submodules/paffy/lib/*.a ${LIBDIR}
	ln -f submodules/paffy/inc/*.h ${INCLDIR}

suball.red:
	cd submodules/red && ${MAKE}
	ln -f submodules/red/bin/Red ${BINDIR}

suball.FASTGA:
	cd submodules/FASTGA && ${MAKE}
	ln -f submodules/FASTGA/FastGA ${BINDIR}
	ln -f submodules/FASTGA/ALNtoPAF ${BINDIR}
	ln -f submodules/FASTGA/FAtoGDB ${BINDIR}
	ln -f submodules/FASTGA/GIXmake ${BINDIR}
	ln -f submodules/FASTGA/GIXrm ${BINDIR}

subclean.%:
	cd submodules/$* && ${MAKE} clean

##
# docker
##
docker: Dockerfile
	-docker rmi -f ${name}:latest
	docker build --network=host -t ${name}:${tag} . --build-arg CACTUS_COMMIT=${git_commit}
	docker tag ${name}:${tag} ${name}:latest
	docker tag ${name}:${tag} evolvertestdocker/cactus:latest

# Requires ~/.dockercfg
push: docker
	docker push ${name}:${tag}
	docker push ${name}:latest

push_only:
	docker push ${name}:${tag}
	docker push ${name}:latest

binRelease:
	./build-tools/makeBinRelease
	CACTUS_LEGACY_ARCH=1 ./build-tools/makeBinRelease

srcRelease:
	./build-tools/makeSrcRelease
