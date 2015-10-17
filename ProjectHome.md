The **SHEL** (SHallow-waters numerical modEL) is a **finite volume**, **free-surface**, **variable bottom**, **shallow-waters equations numerical solver**.

The SHEL is coded in **Matlab** with a built-in graphical interface for loading, editing and saving of simulation parameters and forcings and also for running, visualizing and exporting images (**eps**, **png**) and movies (**avi**).

![http://farm5.static.flickr.com/4119/4910274195_383cdf50b4.jpg](http://farm5.static.flickr.com/4119/4910274195_383cdf50b4.jpg)

The code is compact, efficient and **extensible**, meaning that developers can easily replace the core solver files with custom numerical schemes and can even contribute to the stack of available numerical schemes.

The SHEL, by default, uses an **Arakawa C** grid type over a **land-mask** with a **second-order accurate** in time and space leapfrog and central differences schemes for the **momentum equations** and a **first-order accurate** upwind scheme for the **tracer equation**. **Dirichelet**, **Neumann** and **Sommerfeld** type conditions were implemented as boundary conditions. It is fairly easy to replace these numerical methods with others.

The SHEL is **fully documented** and is **open-source**. The documentation consists of a [user's guide](UsersManual.md), a [technical guide](TechnicalGuide.md), and a [developer's guide](DevelopersManual.md). The software package includes the source-code and is available in the [downloads section](http://code.google.com/p/shel/downloads/list).

If you use SHEL in your work please **cite** the scientific documentation as follows
```
Riflet, G., 2010. SHEL, a Shallow-Water Numerical 
Model: Scientific Documentation. MARETEC, Instituto Superior Técnico, 
Universidade Técnica de Lisboa.
```

To contact the author, please send an email to
`guillaume.riflet at gmail.com`

<a href='http://www.youtube.com/watch?feature=player_embedded&v=90FipZAmCUE' target='_blank'><img src='http://img.youtube.com/vi/90FipZAmCUE/0.jpg' width='425' height=344 /></a>