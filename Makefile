VERSION=$(shell python3 -c "import quadpy; print(quadpy.__version__)")

default:
	@echo "\"make publish\"?"

README.rst: README.md
	cat README.md | sed 's_<img src="\([^"]*\)" width="\([^"]*\)">_![](\1){width="\2"}_g' > /tmp/README.md
	pandoc /tmp/README.md -o README.rst
	python3 setup.py check -r -s || exit 1

# https://packaging.python.org/distributing/#id72
upload: setup.py README.rst
	# Make sure we're on the master branch
	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "master" ]; then exit 1; fi
	rm -f dist/*
	python3 setup.py bdist_wheel --universal
	gpg --detach-sign -a dist/*
	twine upload dist/*

tag:
	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "master" ]; then exit 1; fi
	@echo "Tagging v$(VERSION)..."
	git tag v$(VERSION)
	git push --tags

publish: tag upload

clean:
	rm -f README.rst
