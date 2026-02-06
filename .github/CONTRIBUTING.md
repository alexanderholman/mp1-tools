<!-- Written by Alexander Holman 19/08/2019 -->
<!--
This document is intended to contain an overview or guide on what is expected of anyone that contributes in any fashion 
to the repository. It is written to be used with most internal projects, but as per the other documents feel free to 
adapt it to fit your needs, adding, modifying and removing as necessary.
-->
# Contributing
Before you start contributing to this project please take a moment to read this document. If something is unclear or 
needs clarification, just ask. If you aren't a member of the internal team, this can be done by opening a new issue 
a labelling it as "Type: Question", someone will respond as soon as possible.

## Agreements
By contributing to this project you agree to adhere to the our [CODE OF CONDUCT](/CODE_OF_CONDUCT.md).

## Issues

### Security/Vulnerabilities
Please **DO NOT** create issues relating to security or vulnerabilities via GitHub's issue tracker, instead please do 
so in person or by emailing us directly at [security@holman.org.uk](mailto:security@holman.org.uk).

### Everything else
For all other issues not related to security (including `bugs`, `tasks`, `feature requests`, `questions` and 
`brainstorming`) raise an issue under [GitHub's issue tracker](https://github.com/alexanderholman/{{repo-name}}/issues) 
using the template most fitting.

Before create a new issue, please take a moment to search whether or not a similar or the same issue has already been 
created, and continue to use that issue should it exist. In the event the issue has been closed re-open it.

Please take the time to read through the templated content and be sure to fill in all the information requested. Before 
submitting the issue please review the checklist at the end of the issue. If complete fill the square brackets 
`[ ] - ...` with a captial "X" like so; `[X] - ...`.

#### Types of issues
* **Bug** For issues relating to [bugs](https://en.wikipedia.org/wiki/Software_bug), exceptions, errors or anything 
else where there is something that shouldn't happen, but it's happening! Please take a moment to ensure the issue isn't 
on your end before raising a ticket.
* **Task** `internal only` Slightly more pressing that feature requests, these are for scenarios where something needs 
doing and it's implementation is required, verses a feature request which is a "nice to have" or a suggestion.
* **Feature Request** For suggestions or requests for addition features or changes, generally they should be for 
extensions to current functionality.
* **Question** If you have a question not address in the [SUPPORT](/SUPPORT.md) document or the [Wiki](https://github.com/alexanderholman/{{repo-name}}/wiki) 
(should one exist) please ask it under this label/template. We will respond as soon as we are able to do so.
* **Brainstorm** Issues to hash out ideas and future possibilities, whether with the aim of creating a feature request 
or a task from it or not.
* **Custom** If none of the above "types" fit what is needed, create a custom/blank issue but be sure to label it 
correctly, setting the `type`, `priority`, any `alerts`, `dev` or `skills` required.

#### Assignment
Unless you explicitly know who should be assigned to an issue, please refrain from assigning anyone other than 
yourself, and leave assignment to whomever is undertaking the role of project manager.

## Code Base Contributions
Contributions to the code base is made via pull requests on the GitHub platform.

The primary languages here at [@alexanderholman](https://github.com/alexanderholman) are [PHP](http://www.php.net/), 
and [TypeScript](https://www.typescriptlang.org/) or [JavaScript](https://developer.mozilla.org/bm/docs/Web/JavaScript) 
where TS isn't appropriate.

Regardless of the language used the criteria for contribution remains the same, for what ever language a *standard* is 
defined and must be adhered to at all times. Code *style* should follow a defined approach. Where appropriate, work 
should be *documented*. All vital code should be *tested*. Failure to meet this criteria will result in the 
contribution being rejected and changes requested.

As a rule of thumb, a test-first approach is advised, documenting through code blocks with exceptions to guides, 
tutorials, design diagrams and specifications.

See language specific details below for information relating to *standards*, *style*, *documentation* and *testing*.

### Contribution workflow
The general workflow followed consist of `fork` > `issue` > `branch` > `code` > `test` > `commit` > `pull request` > 
`code review` > `merge` > `version` > `deploy`.

#### Forking
All development should be completed on your own [Fork](https://help.github.com/articles/fork-a-repo/), this is 
equivalent to [cloning](https://git-scm.com/docs/git-clone) a repository on GitHub's platform. Under no circumstances 
should you push code to your own profile/repository that isn't a fork of the [@alexanderholman](https://github.com/alexanderholman) 
project repository.

#### Branching
Each issue being worked on should have it's own branch named after that branch e.g. for an issue titled "Bug in code" 
with an ID of `123`, a branch should be made named `issue-123/bug-in-code`.

With the exception of the deployment branches, no branch should exist without an issue counterpart.

#### Committing
When committing it's easy in an effort to be more efficient to commit you work in a similar method to below:
```bash
git add . && git commit -m "work done etc"
```

The problem with this is that context can quickly become lost. Please follow the commit guidelines set out by [Udacity](https://udacity.com/) 
the online learning platform, which can be found [here](https://udacity.github.io/git-styleguide/).

Commits ideally should made on a per task basis or similar, there is probably no need for a per line walkthrough. 
Think, in 6 months time when I have forgotten to commit to the changelog and need to update it, what information would 
I need to know to rebuild an overview of the changes per release without having to study every `diff` between releases.

#### Pull Requests
Instead of merging or pushing directly upstream, you should push to your forked repository and use GitHub's 
[Pull Request](https://github.com/alexanderholman/{{repo-name}}/pulls) utility as the mechanism to merge. Prior to 
merging, all code should be subject to the review of another developer, even if it was written as a pair.

##### Code Review
To assist in detection of bugs and malware, aid in the development of clean software and readable code, all PR's 
however small should be subject to a code review.

Code review, ideally should be undertake by someone new to the code whom is familiar with the tools used in it's 
development. The code should be scrutinised by the reviewer until it is understood and until all parties are happy. 
Don't be afraid to question where needed or to seek advice from another individual should you need it.

When scrutinising code, consider providing links to articles or anything that may better explain or evidence your point 
of view. Be kind, people make mistakes, including you; ensure your comments and review adheres to our [CODE OF CONDUCT](/CODE_OF_CONDUCT.md) 
where public or not.

##### Merging
Merging to be done as frequently as possibly to keep merging conflicts as few as possible. In the event of a conflict 
consider working side my side with the other developers involved during their resolution.

We aim to keep the entire history tree, as such the use of merging verse of git concatenation tools should be used.

We should only merge upon the following:
* The code has been reviewed by another.
* The code passes all tests provided.
* The code passes inspection of any third-party services, such as [Scrutinizer](https://scrutinizer-ci.com/) and 
[Travis CI](https://travis-ci.org/).

When merging into a deployment branch, ensure no-one else has an issue with it first, e.g. that an update isn't already 
in progress, as it is possible continuous integration and deployment is in place.

### Versioning
Unless required by the standard, [Semantic Versioning](https://semver.org/) following the 2.0.0 standard.

The **TL;DR** is:
* **Major releases** _**major**.minor.patch-meta_ are incompatible with previous versions. Going from `v1.*` to `v2.*` 
compatibility isn't guaranteed and changes will like be required to integrate.
* **Minor releases** _major.**minor**.patch-meta_ are forward compatible but not backwards (while using new features) 
within the same `major` version. Going from `v1.0.*` to `v1.1.*` should require no code changes to update, however 
going from `v1.1.*` down to `v1.0.*` will required either removal of later release feature dependencies or 
implementation.
* **Patch releases** _major.minor.**patch**-meta_ are both forward and backwards compatible within the same `minor` 
version. Going from either `v1.0.0` to `v1.0.1` and `v1.0.1` to `v1.0.0` should require no code changes.
* **Meta release** _major.minor.patch-**meta**_ are informational releases but essential the same as their patch 
equivalents, they may not be adopted e.g. release candidates may be dropped completely.

The above all assumes that contracts and APIs are adhered to.

Not every commit will be versioned, however all deployments should be of versioned code only. Versioning can be done 
through [GitHub's release feature](https://help.github.com/articles/creating-releases/) or [Git tags](https://git-scm.com/book/en/v2/Git-Basics-Tagging).

With the exception of `meta` versions (e.g. `-rc1` or `-dev`), versioning should be maintained for merged code only.

Only released versioned without a metadata extension or where the metadata extension is `stable` should be considered 
a stable release.

### Deployment
<!-- as deployment varies from project to project this section may need updating -->
<!-- EXAMPLE START -->
While the integration process is continuous, the deployments are subject to evaluation. The cycle can be roughly 
represented as follows:

`release alpha` > `merge to development` > `interal developer test` > `release alpha-rc#` > 
`merge to staging` > `internal test` > `release beta` > `merge to staging` > `closed external test` > 
`release beta-rc#` > `merge to staging` > `public external test` > `release stable` > `merge to production`

The act of merging code into the `production`, `staging`, or `development` branches will trigger an update on Jelastic 
based applications within 1 minute. Assuming zero-downtime deployment is enabled, the updates could become visible once 
any of the servers has completed installation. Otherwise, updates will become visible upon installation on all nodes.

As a rule of thumb, `alpha` and `beta` releases will be deployed to `development` and `staging` retrospectively, 
though this isn't a guarantee and it is possible for a release to have it's own testing environment completely. The 
`production` branch is reserved for `stable` releases only and **SHOULD NEVER** contain anything else, be it 
un-versioned or other.
<!--  EXAMPLE END  -->

### Language Specific

#### PHP
For PHP projects, we follow the [PHP Framework Interop Group (PHP-FIP)](https://www.php-fig.org) 
[PHP Standards Recommendations (PSR)](https://www.php-fig.org/psr/) for standards and [Composer](https://getcomposer.org/) 
for dependency management.

Code should be written in an object oriented manner, test first where possible using [PHPUnit](https://phpunit.de/). 
Consider the use of [dependency injection](http://php-di.org/) where possible, as it has been known to aid in testing 
OO applications.

The framework of choice for PHP projects is [Laravel](https://laravel.com/).

Other useful links:
* [Basic Coding Standards PSR-1](https://www.php-fig.org/psr/psr-1/)
* [Coding Style Guidelines PSR-2](https://www.php-fig.org/psr/psr-2/)
* [Autoloading Standard PSR-4](https://www.php-fig.org/psr/psr-4/)

<!-- completion needed
#### TypeScript / JavaScript
{{js-specifics}}

#### Python
{{python-specifics}}

#### C++
{{c++-specifics}}

### Other
{{other-specifics}}
-->