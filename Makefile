#!/usr/bin/make -f

SHELL = /bin/bash

PREFIX                    := /software/ddd/internal/clinical-filter
GIT_URI                   := ssh://git.internal.sanger.ac.uk/repos/git/ddd/clinical-filter.git

TMPDIR                    := $(shell mktemp -d)
SRCDIR                    := $(TMPDIR)/clinical-filter-$(CLINICAL_FILTER_VERSION)

CLINICAL_FILTER_PREFIX    := $(PREFIX)/clinical-filter-$(CLINICAL_FILTER_VERSION)
CLINICAL_FILTER_CONFIGDIR := $(CLINICAL_FILTER_PREFIX)/config

INSTALL                   := /usr/bin/install
CHMOD                     := Du=rwx,Dg=rwx,Do=rx,Fu=rw,Fg=rw,Fo=r




usage:
	@echo "Usage: CLINICAL_FILTER_VERSION=XXX install" >&2

.PHONY: clinical-filter-version

clinical-filter-version:
	@if [ -n "${CLINICAL_FILTER_VERSION}" ]; then echo "Installing clinical-filter-${CLINICAL_FILTER_VERSION}"; else echo "Usage: CLINICAL_FILTER_VERSION=X.X.X install";exit 1; fi

install: clinical-filter-version $(SRCDIR)
	$(MAKE) -C $(SRCDIR) TMPDIR=$(TMPDIR) CLINICAL_FILTER_VERSION=$(CLINICAL_FILTER_VERSION) add-version clean-srcdir-git install-config install-python clean-tmpdir


$(SRCDIR): $(TMPDIR)/clinical-filter-$(CLINICAL_FILTER_VERSION).zip
	cd $(TMPDIR) && unzip $<

$(TMPDIR)/clinical-filter-$(CLINICAL_FILTER_VERSION).zip:
	git archive --format zip --output $(TMPDIR)/clinical-filter-$(CLINICAL_FILTER_VERSION).zip --remote $(GIT_URI) --prefix clinical-filter-$(CLINICAL_FILTER_VERSION)/ $(CLINICAL_FILTER_VERSION)

add-version:
	echo "def version(): return \"$(CLINICAL_FILTER_VERSION)\"" > $(SRCDIR)/src/main/python/clinicalfilter/version.py

clean-srcdir-git:
	find $(SRCDIR)/src/main/ | grep "/\.gitignore"$ | xargs -I '{}' rm {}

install-config: $(SRCDIR)/config/tags.txt $(SRCDIR)/config/filters.txt
	$(INSTALL) -d -m 0755 $(CLINICAL_FILTER_CONFIGDIR)
	$(INSTALL) -m 0644 $^ $(CLINICAL_FILTER_CONFIGDIR) 

install-python:
	rsync -rp --chmod=$(CHMOD) $(SRCDIR)/src/main/python/ $(CLINICAL_FILTER_PREFIX)

clean-tmpdir:
	rm -r $(TMPDIR)

test:
	export PYTHONPATH="src/main/python:src/test/python:${PYTHONPATH}" && python3 -m unittest discover ./src/test/python/clinicalfilter
