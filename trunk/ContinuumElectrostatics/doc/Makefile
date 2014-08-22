PROJECT = ContinuumElectrostatics_manual

default: compile

compile: $(PROJECT).tex refs.bib
	pdflatex   $(PROJECT)
	bibtex     $(PROJECT)
	pdflatex   $(PROJECT)
	pdflatex   $(PROJECT)

clean:
	rm $(PROJECT).aux
