# Linux gcc
COMPILE = gcc -c -fPIC -Wall
LINK = gcc -shared

# Alpha cxx
#COMPILE = cxx -c -O3 -arch host -ev6 -noexception -std gnu
#LINK = cxx -shared

# Alpha gcc
#COMPILE = gcc -c -O3 -fPIC -ev6 -ieee
#LINK = gcc -shared

INCLUDEDEST = ${LIBBINARY}/include/
LIBDEST = ${LIBBINARY}/lib/
LIBNAME = libbinary.so
MAJOR = .1

# Define include and library directories
INCLUDEDIR = -I/usr/include -I/usr/local/include
LIBDIR = -L/usr/local/lib -L/usr/lib -L/usr/X11R6/lib

OBJS = stream.o binary.o visibility.o misc.o idl_binary.o

# Define libraries to use for linking
LIBS = # Don't need any for a shared library

$(LIBNAME)$(MAJOR) : $(OBJS)
	$(LINK) -Wl,-soname,$(LIBNAME)$(MAJOR) $(LIBDIR) -o $(LIBNAME)$(MAJOR) $(OBJS) $(LIBS)

stream.o : stream.c
	$(COMPILE) $(INCLUDEDIR) stream.c

binary.o : binary.c
	$(COMPILE) $(INCLUDEDIR) binary.c

visibility.o : visibility.c
	$(COMPILE) $(INCLUDEDIR) visibility.c

misc.o : misc.c
	$(COMPILE) $(INCLUDEDIR) misc.c

idl_binary.o : idl_binary.c
	$(COMPILE) $(INCLUDEDIR) idl_binary.c

install : $(LIBNAME)$(MAJOR)
	mkdir -p $(INCLUDEDEST)
	mkdir -p $(LIBDEST)
	cp *.h $(INCLUDEDEST); echo "Installed include files in $(INCLUDEDEST)"; fi
	cp $(LIBNAME)$(MAJOR) $(LIBDEST)/;  echo "Installed $(LIBNAME) files in $(LIBDEST)";
	rm $(LIBDEST)/$(LIBNAME) || echo "No existing  $(LIBDEST)/$(LIBNAME)"
	ln -s $(LIBNAME)$(MAJOR) $(LIBDEST)/$(LIBNAME); fi
	echo $(LIBDEST) > /etc/ld.so.conf.d/libbinary.conf
	/sbin/ldconfig

clean:
	rm -f *.o
	rm -f *~
