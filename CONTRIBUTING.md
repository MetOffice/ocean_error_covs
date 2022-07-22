# How to contribute

## Create a Branch 

* Clone the repository 

```
git clone https://github.com/MetOffice/ocean_error_covs.git
```

* Create your local branch 

```
git checkout -b <BranchName>
```

The **BranchName** should have the following format ```<institution>/<github_username>/<issue_type>/<branch_title>```.

* Institution: the institution where you work (e.g. MetOffice)
* Github Username: your username
* Issue Type:
  * feature: for new developments
  * enhancement: for improving features
  * bug: for fixing bugs
* Branch title: what briefly describes your development

An example is provided below: 

```
git checkout -b MetOffice/mo-dcarneir/feature/write_code_guidelines
```

* Push this branch back to repository to store it remotely

```
git fetch origin
git push origin <BranchName>
```

## Merge Requests

All contributions to ```ocean_error_covs``` are made via merges with the master branch. Once the branch is ready to be merged, a pull request should be created in the "Pull requests" section of ```ocean_error_covs```. At least one developer must be assigned to review the pull request. The developer who reviews the pull request is responsbile for checking the changes before approving the merge request.

:warning: Before creating the merge pull request you must run the test cases and ensure that they all match their respective Known Good Outputs (KGOs). This is mandatory and needs to be informed in the merge pull request! 

After merging the branch, do not forget to close your issue related to that merge if you have created one. 

## Code Contributors 

New contributors are encouraged to add their details here as part of their first request. 

* Davi Mignac Carneiro (Met Office, UK) - davi.carneiro@metoffice.gov.uk
* Matthew Martin (Met Office, UK) - matthew.martin@metoffice.gov.uk 
* James While (Met Office, UK) - james.while@metoffice.gov.uk
