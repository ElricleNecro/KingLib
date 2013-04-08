CC=gcc

CFLAG=-std=c99 -g3 -W -Wall -Wshadow -Wcast-qual \
      -Wcast-align -Wsign-compare -Wstrict-prototypes \
      -Wredundant-decls \
      -Wnested-externs \
      -ffloat-store -Wunreachable-code -Wwrite-strings -fPIC
#-DUSE_FENV
INC=-I include/
LFLAG=-fPIC
LINK=
#-lm

TAR=tar
TAR_OPT=-cvzf
TAR_EXT=tar.gz

SRCDIR=src/
OBJDIR=build/
LIBDIR=build/

PREFIX=$$HOME/.local
PYPREFIX=$$HOME/.local/lib/python3.3/site-packages/

PROJ=king
LIBNAME=lib$(PROJ)

SRC=king.c mod.c rand.c rk4.c utils.c
OBJ=$(SRC:.c=.o)
HEA=include/king/cte_phys.h include/king/isotherm.h include/king/king.h include/king/mod.h include/king/rand.h include/king/rk4.h include/king/utils.h

WRAP_SRC=$(PROJ)PYTHON_wrap.c
SWIG_SRC=$(PROJ).i
SWIG_OPT=-Wall -python -py3
WRAP_LIB=_$(PROJ).so

#!all:
#!	Crée la librairie et le module python 3.
#!
all:$(LIBDIR)/$(LIBNAME).so PyWrap

#!PyWrap:
#!	Crée juste le module (a pour dépendance la librairie, équivaut donc à all).
#!
PyWrap:$(LIBDIR)/$(LIBNAME).so $(LIBDIR)/$(WRAP_LIB)

#!lib:
#!	Crée juste la librairie sans le module python.
#!
lib:$(LIBDIR)/$(LIBNAME).so

#!install-all:
#!	Installe la librairie et le module python.
#!	Pensez à jouer avec les variables PREFIX et PYPREFIX pour changer les répertoires d'installation.
#!	Par défaut :
#!		PREFIX=$HOME/.local
#!		PYPREFIX=$HOME/.local/lib/python3.2/site-packages/
#!
install-all:install install-py

#!install:
#!	Installe la librairie en utilisant la variable PREFIX
#!
install:$(LIBDIR)/$(LIBNAME).so
	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/include
	cp $< $(PREFIX)/lib/.
	cp -r include/$(PROJ) $(PREFIX)/include/.
	TMP=$(PREFIX) sed -e "s:HOME:$$TMP:g" king.pc > $(PREFIX)/lib/pkgconfig/king.pc

#!install-py:
#!	Installe la librairie en utilisant la variable PYPREFIX
#!
install-py:$(LIBDIR)/$(WRAP_LIB)
	mkdir -p $(PYPREFIX)/$(PROJ)
	cp $< $(PYPREFIX)/$(PROJ)/.
	cp $(LIBDIR)/$(SWIG_SRC:.i=.py) $(PYPREFIX)/$(PROJ)/.

#$(OBJDIR):
#	@mkdir -p $@
#
#$(LIBDIR):
#	@mkdir -p $@

#$(OBJDIR)
$(OBJDIR)/%.o:$(SRCDIR)/%.c $(HEA)
	$(CC) $(CFLAG) $(INC) -c $< -o $@

#$(LIBDIR)
$(LIBDIR)/$(LIBNAME).so:$(foreach x, $(OBJ), $(OBJDIR)/$(x)) $(HEA)
	$(CC) $(LFLAG) -shared $(foreach x, $(OBJ), $(OBJDIR)/$(x)) -o $@ $(LINK)

$(WRAP_SRC):$(SWIG_SRC) $(LIBDIR)/$(LIBNAME).so
	swig -Wall -python -py3 -I/usr/include -outdir $(LIBDIR) -o $@ $<

$(LIBDIR)/$(WRAP_LIB):$(WRAP_SRC)
	$(CC) $(CFLAG) `pkg-config --cflags python3` $(INC) -I../ -fPIC -shared -o $@ $< -L $(LIBDIR) -l$(PROJ) $(LINK) `pkg-config --libs python3`

#!help:
#!	Affiche cette aide.
#!
help:
	@grep "^#!" Makefile | sed -e 's/#!//'

#!test:
#!	Éxecute un test afin de vérifier le bon fonctionnement du module python.
#!
test:
	./test_king-py3.py
	@echo "Utilisez votre visualiseur d'image favori pour regarder l'image test_py3KING.png"

#!doc:
#!	Crée la documentation de la librairie (suffisante pour le module python).
#!
doc:
	doxygen -g Doxyfile

#!tar:
#!	Crée une archive des élèments essentiel à la compilation et au test du code.
#!
tar:
	$(TAR) $(TAR_OPT) $(PROJ).$(TAR_EXT) ../KingLib/src/* ../KingLib/include/$(PROJ)/* ../KingLib/Makefile ../KingLib/test_king-py3.py ../KingLib/test_king-py2.py ../KingLib/King_dim.log ../KingLib/TODO ../KingLib/$(SWIG_SRC) ../KingLib/INSTALL ../KingLib/config.in ../KingLib/Doxyfile ../KingLib/logo.png

tar-all:
	$(TAR) $(TAR_OPT) Perso_$(PROJ).$(TAR_EXT) ../KingLib/*

#!uninstall:
#!	Désinstalle les élèments de la librairie et du module python en usant des variables PREFIX et PYPREFIX.
#!	Attention, si vous les avez modifiez, pensez à les ajuster.
#!
uninstall:
	$(RM) -r $(PYPREFIX)/$(PROJ)
	$(RM) -r $(PREFIX)/include/$(PROJ)
	$(RM) $(PREFIX)/lib/$(LIBDIR)/$(LIBNAME).so

#!clean:
#!	Nettoyage simple. Efface les fichiers objets.
#!
clean:
	$(RM) $(foreach x, $(OBJ), $(OBJDIR)/$(x))

#!clean-all:
#!	Nettoyage complet. Efface les fichiers objets, les fichiers associés au module python, et les librairies.
#!
clean-all:clean
	$(RM) $(LIBDIR)/$(LIBNAME).so $(LIBDIR)/$(WRAP_LIB) $(WRAP_SRC) $(LIBDIR)/$(SWIG_SRC:.i=.py)

#!clean-doc:
#!	Nettoyage spécifique. Efface la documentation.
#!
clean-doc:
	$(RM) -r Doc/

#!clean-dist:
#!	Super-nettoyage. Efface tout en appelant les cibles clean précédentes.
#!
clean-dist:clean-all clean-doc

