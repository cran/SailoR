# Version 1.1

* Improved visibility of the Sailor diagrams by plotting points on top of ellipses. The kind of character used depends on the type of plot produced (centered or uncentered) in order to improve visibility, too.

* Added new diagnostics (eccentricity of the ellipses and congruence coefficient of the leading eigenvector) to improve the interpretability of results, particularly the rotation angle of EOFs (congruence coefficient) and its sensitivity to degeneragy of eigenvalues (degeneracy would be higher with eccentricities closer to zero).

* A new synthetic dataset and new examples were added. The main objective of adding this synthetic dataset is to show the response of the Sailor diagram to different kinds of errors.

# Version 1.1.1

* A new parameter was added to SailoR.plot() function: line type. The user can now choose the type of line to be used to plot the diagram.

# Version 1.1.2

* A new parameter was added to SailoR.Plot() function: bias symbol type. The user can now choose the type of symbol to be used to plot the bias.

* The Ensembles dataset and its manual page were updated: A new model was added to it. 

* A new example was added to SailoR.Table(). This new example shows how to generate the table in LaTeX.

# Version 1.2

* All users should update their SailoR package to this new version. A bug affecting ONLY the estimations of the RMSE values has been solved. The rest of the terms in the indices and output from the UVError function were correct. Thanks to Kwun Yip (Samuel) Fung for reporting us about it. 

