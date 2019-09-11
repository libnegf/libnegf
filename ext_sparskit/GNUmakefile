.SUFFIXES: 

vpath % $(SRCDIR)

ifeq ($(strip $(OBJDIR)),)
  OBJDIR:="."	
endif

AR = ar -rc 
TARGET = $(OBJDIR)/libsparskit.a

SOURCES = \
    BLASSM/blassm.f \
    BLASSM/zblassm.f \
    BLASSM/matvec.f \
    BLASSM/zmatvec.f \
    FORMATS/formats.f \
    FORMATS/zformats.f \
    FORMATS/unary.f \
    FORMATS/zunary.f \
    FORMATS/sort.f \
    INFO/infofun.f \
    INFO/zinfofun.f \
    INOUT/inout.f \
    INOUT/zinout.f \
    ORDERINGS/ccncommon.f \
    ORDERINGS/ccn.f \
    ORDERINGS/zccn.f \
    ORDERINGS/zcolor.f \
    ORDERINGS/zdsepart.f
OBJS = $(SOURCES:.f=.o)	

DIRS := $(sort $(dir $(OBJS)))
DIR_TIMESTAMPS = $(addsuffix _TIMESTAMP_,$(DIRS))

.PHONY: all

all: $(TARGET) 

$(DIR_TIMESTAMPS):
	mkdir -p $(@D)
	touch $@

$(TARGET): $(OBJS) 
	$(AR) $(TARGET) $(OBJS) 

%.o: %.f $(DIR_TIMESTAMPS)
	$(FXX) $(FXXOPT) -o $@ -c $<


.PHONY: clean
clean:
	rm -rf $(DIRS)

.PHONY: distclean
distclean: clean
	rm -f $(TARGET)
