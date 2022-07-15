
# Set ROOT environment variables
ROOTSYS = /afs/rhic.bnl.gov/star/ROOT/5.34.30/.sl64_gcc482/rootdeb
ROOTCFLAGS = -pthread -std=c++11 -Wno-deprecated-declarations -m32 -I/afs/rhic.bnl.gov/star/ROOT/5.34.30/.sl64_gcc482/rootdeb/include
ROOTLIBS = -L/afs/rhic.bnl.gov/star/ROOT/5.34.30/.sl64_gcc482/rootdeb/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic
ROOTGLIBS = -L/afs/rhic.bnl.gov/star/ROOT/5.34.30/.sl64_gcc482/rootdeb/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic
#<
