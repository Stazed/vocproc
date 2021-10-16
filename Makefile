
# By default plugin will be installed in /usr/local/lib/lv2.
# You can change that by changing value of INSTALL_DIR

DESTDIR=
INSTALL_DIR = $(DESTDIR)/usr/lib/lv2


# By default (because of the built in vocoder functionality) VocProc has two inputs and
# one output. If you're using VocProc in Ardour, that means that you need to apply it 
# always to stereo track in whick one channel will be voice input and other carrier (empty
# if vocoder is not used)
#
# Setting DEFINES to -DNO_VOCODER removes built in vocoder so the final plugin has only
# one input and one output. Use NO_VOCODER if you plan to use it mostly in Ardour.If you use 
# it in hosts like Ingen where you can wire it by your wishes, you can leave the default setting.

#CPPFLAGS=$(shell pkg-config --cflags gtkmm-3.0)
CPPFLAGS=$(shell pkg-config --cflags --libs glibmm-2.4 giomm-2.4 gtkmm-3.0)
# do not edit after this point

BUNDLE = vocproc.lv2

all: $(BUNDLE)

$(BUNDLE): manifest.ttl vocproc.ttl vocproc.so vocproc_gui.so vocproc_gui.ui 
	rm -rf $(BUNDLE)
	mkdir $(BUNDLE)
	cp $^ $(BUNDLE)

vocproc.so: vocproc.cpp vocproc.peg
	g++ $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) -shared -fPIC -DPIC $(DEFINES) vocproc.cpp `pkg-config --cflags --libs lv2-plugin fftw3` -g -lm -o vocproc.so

vocproc_gui.so: vocproc_gui.cpp vocproc.peg
	g++ $(LDFLAGS) $(CPPFLAGS) $(CFLAGS) -shared -fPIC -DPIC $(DEFINES) vocproc_gui.cpp `pkg-config --cflags --libs lv2-gui` -o vocproc_gui.so

vocproc.peg: vocproc.ttl
	lv2peg vocproc.ttl vocproc.peg

install: $(BUNDLE)
	mkdir -p $(INSTALL_DIR)
	rm -rf $(INSTALL_DIR)/$(BUNDLE)
	cp -R $(BUNDLE) $(INSTALL_DIR)

uninstall:
	rm -rf $(INSTALL_DIR)/$(BUNDLE)

clean:
	rm -rf $(BUNDLE) vocproc.so vocproc_gui.so vocproc.peg

