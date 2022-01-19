# How to contribute

## Raise Issues

Report bugs and request enhancements by raising an issue in the ```ocean_error_covs``` issues section. If reporting an issue, try to include information on how to reproduce the bug/error. It is also important to define the assignee(s) and label the issue, indicating whether it is related to a code enhancement, bug or question/help. 

## Create a Branch 

If the issue involves the creation of a branch, please follow the protocol below after raising the issue: 

* Clone the repository 

```
git clone https://github.com/MetOffice/ocean_error_covs.git
```

* Create your local branch 

```
git checkout -b <BranchName>
```

The **BranchName** should have the following format ```<github_username>_issue#<issue_number>_<feature>```. An example is provided below: 

```
git checkout -b mo-dcarneir_issue#1_code_guidelines
```

* Push this branch back to repository to store it remotely

```
git fetch origin
git push origin <BranchName>
```

After following the steps above you are ready to make changes to your branch.

## Merge Requests

All contributions to ```ocean_error_covs``` are made via merges with the master branch. Once the branch is ready to be merged, a pull request should be created in the "Pull requests" section of ```ocean_error_covs```. At least one developer must be assigned to review the pull request. The developer who reviews the pull request is responsbile for checking the changes before approving the merge request. After merging the branch, do not forget to close your issue related to that merge.  

New contributors are encouraged to add their details to the code contributors section of this document as part of their first request. 

## Code Contributors 

* Davi Mignac Carneiro (Met Office, UK) - davi.carneiro@metoffice.gov.uk
* Matthew Martin (Met Office, UK) - matthew.martin@metoffice.gov.uk 
* James While (Met Office, UK) - james.while@metoffice.gov.uk
