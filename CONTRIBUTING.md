# Contributing guidelines

We thank you in advance :thumbsup: :tada: for taking the time to contribute, whether with *code* or with *ideas*, to the project.


## Did you find a bug?

* Ensure that the bug was not already reported by [searching under Issues].

* If you're unable to find an (open) issue addressing the problem, [open a new one]. Be sure to prefix the issue title with **[BUG]** and to include:

  - a *clear* description,
  - as much relevant information as possible, and
  - a *code sample* or an (executable) *test case* demonstrating the expected behaviour that is not occurring.

## How to work on a new feature/bug

Create an issue on Github or you can alternatively pick one already created.

Assign yourself to that issue.

Discussions on how to proceed about that issue take place in the comment section on that issue.

Some of the work might have been done already by somebody, hence we avoid unnecessary work duplication and a waste of time and effort. Other reason for discussing the issue beforehand is to communicate with the team the changes as some of the features might impact different components, and we can plan accordingly.

## How we work with Git

All work should take place in a dedicated branch with a short descriptive name.

Use comments in your code, choose variable and function names that clearly show what you intend to implement.

Once the feature is done you can request it to be merged back into `master` by making a Pull Request.

Before making the pull request it is a good idea to rebase your branch to `master` to ensure that eventual conflicts with the `master` branch is solved before the PR is reviewed and we can therefore have a clean merge.


### General stuff about git and commit messages

In general it is better to commit often. Small commits are easier to roll back and also makes the code easier to review.

Write helpful commit messages that describes the changes and possibly why they were necessary.

Each commit should contain changes that are functionally connected and/or related. If you for example want to write _and_ in the first line of the commit message this is an indicator that it should have been two commits.

Learn how to select chunks of changed files to do multiple separate commits of unrelated things. This can be done with either `git add -p` or `git commit -p`.


#### Helpful commit messages

The commit messages may be seen as meta-comments on the code that are incredibly helpful for anyone who wants to know how this piece of software is working, including colleagues (current and future) and external users.

Some tips about writing helpful commit messages:

 1. Separate subject (the first line of the message) from body with a blank line.
 2. Limit the subject line to 50 characters.
 3. Capitalize the subject line.
 4. Do not end the subject line with a period.
 5. Use the imperative mood in the subject line.
 6. Wrap the body at 72 characters.
 7. Use the body to explain what and why vs. how.

For an in-depth explanation of the above points, please see [How to Write a Git Commit Message](http://chris.beams.io/posts/git-commit/).


### How we do code reviews

A code review is initiated when someone has made a Pull Request in the appropriate repo on github.

Work should not continue on the branch _unless_ it is a [Draft Pull Request](https://github.blog/2019-02-14-introducing-draft-pull-requests/). Once the PR is marked ready the review can start.

The initiator of the PR should recruit a reviewer that get assigned reviewer duty on the branch. 

Other people may also look at and review the code.

A reviewers job is to:

  * Write polite and friendly comments - remember that it can be tough to have other people critizising your work, a little kindness goes a long way. This does not mean we should not comment on things that need to be changed of course.
  * Read the code and make sure it is understandable
  * Make sure that commit messages and commits are structured so that it is possible to understand why certain changes were made.

Once the review is positive the Pull Request can be _merged_ into `master` and the feature branch deleted.


----

Thanks again.

[searching under Issues]: https://github.com/NBISweden/GAAS/issues?utf8=%E2%9C%93&q=is%3Aissue%20label%3Abug%20%5BBUG%5D%20in%3Atitle
[open a new one]: https://github.com/NBISweden/GAAS/issues/new?title=%5BBUG%5D
