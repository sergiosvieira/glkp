TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -L./glpk-4.60/lib -lglpk
INCLUDEPATH += . ./glpk-4.60/include
QMAKE_CXXFLAGS += -nologo -W3 -O2 -Zi
SOURCES += main.cpp
