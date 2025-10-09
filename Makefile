DOCKER_IMAGE := ydmt/dreem
VERSION := $(shell git describe --always --dirty --long)
PYPI_PASSWORD := $(shell cat ~/.pypi_pass.txt)

default: 
	python setup.py install

pytest:
	pytest test -v

init:
	pip install -r requirements.txt

build-image:
	docker build .
		-f ./Dockerfile
		-t $(DOCKER_IMAGE):$(VERSION)

push-image:
	docker push $(DOCKER_IMAGE):$(VERSION)

upgrade-dependencies:
	pip uninstall -y seismic-graph
	rm -f requirements.txt
	pip freeze > requirements.txt
	python setup.py install

push_to_pypi:
	rm -fr dist
	python3 -m build
	twine upload -r pypi dist/* --user yvesmartindestaillades --password $(PYPI_PASSWORD)
	rm -fr dist

docs_update:
	rm -rf docs/_*
	rm -rf docs/.doctrees
	rm -rf docs/*.html
	rm -rf docs/*.inv
	rm -rf docs/*.js
	sphinx-build docs/source docs/build docs/* -b html

docs_clear:
	rm -rf docs/_*
	rm -rf docs/.doctrees
	rm -rf docs/*.html
	rm -rf docs/*.inv
	rm -rf docs/*.js

# Visual regression testing targets
visual-baseline:
	python test/visual_regression/test_runner.py baseline

visual-test:
	python test/visual_regression/test_runner.py current

visual-compare:
	python test/visual_regression/compare.py list

visual-compare-open:
	python test/visual_regression/compare.py open-all

visual-compare-diff:
	python test/visual_regression/compare.py diff

visual-approve:
	python test/visual_regression/compare.py approve

visual-help:
	@echo "Visual Regression Testing Commands:"
	@echo "  make visual-baseline       Generate baseline plots (before code changes)"
	@echo "  make visual-test          Generate current plots (after code changes)"
	@echo "  make visual-compare       List all test cases"
	@echo "  make visual-compare-open  Open all tests in browser for visual comparison"
	@echo "  make visual-compare-diff  Show JSON structure differences"
	@echo "  make visual-approve       Approve current as new baseline (after review)"
	@echo ""
	@echo "Typical workflow:"
	@echo "  1. make visual-baseline     # Before making changes"
	@echo "  2. [Make code changes]"
	@echo "  3. make visual-test         # After making changes"
	@echo "  4. make visual-compare-open # Review changes in browser"
	@echo "  5. make visual-approve      # If changes look good"
