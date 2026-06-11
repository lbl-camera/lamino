.PHONY: help clean html pdf
.DEFAULT_GOAL := help

help:
	@echo "Tomocam Documentation Build"
	@echo ""
	@echo "Available targets:"
	@echo "  html       - Build HTML documentation"
	@echo "  pdf        - Build PDF documentation"
	@echo "  clean      - Remove build artifacts"
	@echo "  doxygen    - Generate Doxygen XML"
	@echo "  serve      - Start local documentation server"
	@echo ""

clean:
	rm -rf _build/ doxygen/

doxygen:
	doxygen Doxyfile

html: doxygen
	sphinx-build -b html . _build/html

pdf:
	sphinx-build -b latex . _build/latex
	cd _build/latex && make

serve: html
	@echo "Starting documentation server at http://localhost:8000"
	cd _build/html && python -m http.server 8000

install-deps:
	pip install -r requirements.txt
