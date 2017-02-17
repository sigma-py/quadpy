VERSION=$(shell python -c "import quadpy; print(quadpy.__version__)")

# Make sure we're on the master branch
ifneq "$(shell git rev-parse --abbrev-ref HEAD)" "master"
$(error Not on master branch)
endif

default:
	@echo "\"make publish\"?"

README.rst: README.md
	cat README.md | sed 's_<img src="\([^"]*\)" width="\([^"]*\)">_![](\1){width="\2"}_g' > /tmp/README.md
	pandoc /tmp/README.md -o README.rst
	python setup.py check -r -s || exit 1

upload: setup.py README.rst
	python setup.py sdist upload --sign

tag:
	@echo "Tagging v$(VERSION)..."
	git tag v$(VERSION)
	git push --tags

publish: tag upload

clean:
	rm -f README.rst
