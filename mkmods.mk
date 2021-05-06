# Ensure the $(MODDIR) exists
$(shell test -d $(MODDIR) || mkdir -p $(MODDIR))

SUFFIXES = .$(FC_MODEXT) _mod.$(FC_MODEXT)
.F90.$(FC_MODEXT) .F90_mod.$(FC_MODEXT) .f90.$(FC_MODEXT) .f90_mod.$(FC_MODEXT):
	$(PPFCCOMPILE) -c $<
	@cp $(MODDIR)/$@ .

CLEANFILES = *.$(FC_MODEXT) $(BUILT_SOURCES:%=$(MODDIR)/%) *__genmod.$(FC_MODEXT) *__genmod.f90
