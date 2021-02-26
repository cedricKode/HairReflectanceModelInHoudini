QT += core gui opengl

TARGET = cyHairConvert
OBJECTS_DIR=obj

isEqual(QT_MAJOR_VERSION, 5) {
    cache()
    DEFINES +=QT5BUILD
}

MOC_DIR=moc
CONFIG-=app_bundle

INCLUDEPATH += \
    src \
    include \
    models

DESTDIR=./

SOURCES += \
    src/Convert.cpp \
    src/main.cpp


HEADERS += \
    include/Convert.h \
    include/cyHairFile.h

OTHER_FILES += \
    models/*
