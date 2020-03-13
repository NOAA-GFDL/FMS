# Contributing To libFMS

Thank you for taking time to contribute.

libFMS is a software framework for supporting the efficient development,
construction, execution, and scientific interpretation of atmospheric, oceanic,
and climate system models.

Contributions to this project are released to the public under the
[projects open source license](LICENSE.md).

What follows is a set of guidelines and how-tos for contributing to libFMS.
These are guidelines, not rules.  Use your best judgement and feel free to
propose changes to this document in a pull request.


Table of Contents
* [Code of Conduct](#code-of-conduct)
* [Quick Start Workflow](#quick-start-workflow)
* [Support](#support)
* [Issues](#issues)
* [Pull Requests](#pull-requests)
* [Tests](#tests)
* [Styleguides](#styleguides)
* [Maintainer Help](#maintainer-help)

## Code of Conduct

This project and everyone participating in it is governed by the
[Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to
uphold this code. Please report unacceptable behavior to
[gfdl.climate.model.info@noaa.gov](mailto:gfdl.climate.model.info@noaa.gov).

## Quick Start Workflow

1. Create an issue that describes the change, or feature request
2. Fork the project
3. Create a feature branch off an up-do-date `master` branch
4. Update the tests and code
5. Push the commits to your fork
6. Submit a pull request to the `master` branch

## Support

To get support, open a GitHub issue using the support issue template.  Please be
as descriptive as possible.  See the [Issues](#issues) section for more details.

Support for this project is primarily accomplished via this project’s community.
Additional support may be offered from related communities
(e.g. [MOM6](https://github.com/NOAA-GFDL/MOM6) or by members of the GFDL
[Modeling Systems](https://www.gfdl.noaa.gov/modeling-systems) Group.  The
members of the Modeling Systems group are the main developers and maintainers of
this project.  Modeling Systems is a small group, and our ability to offer
support is limited.  Please be patient when requesting support.

## Issues

When contributing to this project, please first open an issue using the template
that best matches the change type.

| Template        | Description                                                        |
| --------------- | ------------------------------------------------------------------ |
| Bug             | Report a bug, or unexpected behavior                               |
| Feature Request | A requested feature or change to the code                          |
| Support         |  Request for support.  Please see the [Support](#support) section. |

The issue title should be short and descriptive.  The issue description should
be clear and concise.  Include enough information to help others reproduce the
issue, or understand the change requested.  Use
[GitHub flavored markdown](https://guides.github.com/features/mastering-markdown/)
to format the message, and include links to references.

## Pull Requests

Submit pull requests for bug fixes, improvements, including tests, or alerts to
particular problems in the code.  We perform merge requests based on an internal
GFDL schedule that addresses the needs of the GFDL scientists.  This release
schedule is better described in the Release Schedule sections.  We follow the
[GitHub fork-pull request workflow](https://guides.github.com/activities/forking/)
workflow, briefly described in the [Quick Start Workflow](#quick-start-workflow)
section.  When creating the pull request, please select a template that best
matches the type of changes.

Please keep the changes in a single pull to be as small as possible to help
reviewer(s) quickly evaluate changes.  If you do have a large update, try to
split the update into small, logical pull requests.  If the update includes code
refactoring, please submit a pull request that includes just the code refactoring.

Once a pull request is created, a member of the Modeling Systems group will
review the changes, and, if necessary, will work with the author of the pull
request to modify their code changes. Note that merging pull requests is
contingent on several factors, including the ability to support the code changes
long-term, portability, and the scope of the impact on the code base. Therefore,
Modeling Systems does not guarantee that all pull requests will be accepted,
even if the changes pass the initial testing phases, and are otherwise correct.

## Tests

FMS uses TravisCI and gitlab-CI to run build tests for libFMS.  Users may create
unit tests, code coverage tests, and regression tests for new and existing code
in yaml (.yml) files.  Github provides a guide
(https://help.github.com/en/articles/about-continuous-integration) for
implementing CI tests.

## Code Style
Code updates should follow the coding style for the project, contained in
[CODE_STYLE.md](CODE_STYLE.md).

## Release Schedule

Releases will be tagged using the format `yyyy.rr[.pp]`, where `yyyy` is the
4-digit year, `rr` is the 2-digit release number, and `pp` is the 2-digit patch
number.  Preliminary releases mean for testing (i.e., code that is still under
development) will be marked `yyyy.rr.alpha.pp` or `yyyy.rr.beta.pp`. Alpha tags
mark code updates that are intended for developers to include in their baseline
regression tests to determine whether the code contains bugs not encountered
during baseline testing. Beta tags are intended for a wider audience of
developers and end users to ensure that their simulations run as expected and
reproduce prior results.
