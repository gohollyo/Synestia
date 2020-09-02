# Abstract
A synestia is a newly discovered planetary object that can evolve into a planet-moon system. If the Moon-forming giant impact resulted in a synestia, then the Moon would have formed within Earth. Moons condense and grow inside a synestia and a synestia eventually cools and transitions into a planet. However interesting this implication may be, much of the public and many researchers are not familiar with synestias and also have diffculty understanding the complex physics which dictate a synestia's behavior. To make learning about synestias more accessible and intuitive, I have developed a set of online Jupyter Notebooks (text, media, and code) for a user to visually explore various physics-based aspects of synestias. In this set of notebooks, a user will learn about: what synestias conceptually look like, the gravitational and extreme thermal profiles of a synestia, the criteria for an object to be a synestia, how synestias form, the balance of forces within a synestia, how gas drag affects moon formation within synestias, and how angular momentum affects the orbital dynamics of moons forming within synestias.

# Introduction
Research is difficult; explaining it is even more so. How can researchers distill complex, nuanced research results into a friendly and easy-to-understand format for others? <i>Jupyter Notebooks</i>, documents that are part text, media, and code that can be placed online ([Kluyver et al., 2016](http://ebooks.iospress.nl/publication/42900)), are an effective visualization educational tool for conveying ideas and research results ([Barba & Forsyth, 2018](https://jose.theoj.org/papers/10.21105/jose.00021.pdf); [Freeman et al., 2014](https://www.pnas.org/content/111/23/8410?lipi=urn%3Ali%3Apage%3Ad_flagship3_pulse_read%3B5ujlJ92ZQgC6PXO%2BbkuCcQ%3D%3D&utm_source=SwitchUp&utm_medium=Blog)). Part of what makes Jupyter Notebooks effective learning tools is their interactive python widgets, <i>ipywidgets</i>. These widgets can be manipulated so that a user can interact with the Jupyter Notebook. Widgets are commonly used in this set of notebooks to allow a user to manipulate plots. A user can adjust the inputs for an analysis, and changes in the plots are shown in real time. This allows the user to explore concepts visually and gain an intuitive understanding of how various parameters affect a particular outcome. Jupyter Notebooks also have a higher readability and approachability than scientific papers, which tend to be less inviting for non-specialists.

## Purpose
The goal of this particular set of Jupyter Notebooks is to help people at the undergraduate level and above understand a newly discovered type of planetary object called a <i>synestia</i>. This discovery is not well understood by those outside of the lunar research community, although it has received a decent amount of publicity [[TED](https://www.ted.com/talks/sarah_t_stewart_where_did_the_moon_come_from_a_new_theory?language=en) (Stewart, 2019); [Scientific American](https://www.ted.com/talks/sarah_t_stewart_where_did_the_moon_come_from_a_new_theory?language=en) (Lock & Stewart, 2019)]. Every new idea takes a while to absorb into the public sphere, still, those who do know about synestias struggle with grasping the physics of their internal dynamics. This is to be expected; we have no real-life reference for synestias (e.g. telescopic images). We are not familiar with this new type of planetary object, so we don't have any intuition for how a synestia might behave. The physics that can explain how synestias work are unfamiliar to us; what we observe on Earth on a daily basis is very different from how those same objects would behave within a synestia's oblate, high-pressure, high-temperature environment. These notebooks seek to bridge that gap -- to help you analyze and gain intuition about the inner workings of synestias so that they seem less esoteric.

## REBOUND
Simulations in this paper made use of the REBOUND code which is freely available at [GitHub](http://github.com/hannorein/rebound) ([Rein & Liu, 2012](https://www.aanda.org/articles/aa/abs/2012/01/aa18085-11/aa18085-11.html)). REBOUND is an open-source software that computes force interactions (gravity in particular) among spherical particles (can represent stars, planets, moons, ring or dust particles). It is efficient and easy to customize. The integrator in the simulations is IAS15, which uses an adaptive time step ([Rein & Spiegel, 2015](https://academic.oup.com/mnras/article/446/2/1424/2892331)).

## Setting up These Jupyter Notebooks on Your Computer
To launch these notebooks and run them in a web browser via Binder, click [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gohollyo/Synestia/master). Note that it may take a few minutes to load.

In order to use these Jupyter Notebooks on your own computer, you will need Python, Jupyter, and REBOUND.

The easiest way to download both Jupyter and Python in one go is by downloading the [Anaconda distribution](https://www.anaconda.com/products/individual). Go to the link, scroll to the bottom of the webpage, and click on an installer according to your computer's operating system. Follow the installation instructions provided on the download webpage.

To download the Python version of REBOUND, in your terminal type either:

`pip install rebound`

or

`conda install -c conda-forge rebound`

Once you have downloaded Python, Jupyter, REBOUND, and these notebooks onto your computer, you can run Jupyter by clicking on the Jupyter icon in the Anaconda app or by typing in the terminal:

`jupyter notebook`

## What's to Come
In the Jupyter notebooks, you are guided through the various dynamics and attributes of synestias, mainly via interactive figures.  
1. <b>What are Synestias? How Were Synestias Discovered?</b> (Giant Impact Origin of the Moon)
2. <b>Physical Properties of Synestias</b> (Temperature, Pressure, Density, Gravity)
3. <b>What Determines Whether a Planetary Structure is a Synestia?</b> (Corotation Limit)
4. <b>How Do Synestias Form?</b> (Giant Impacts)
5. <b>What are the Forces Acting Within Synestias?</b> (Hydrostatic Equilibrium)
6. <b>How Does Gas in a Synestia Affect Moon Formation?</b> (Gas Drag)
7. <b>How Do Moonlets Forming Within Synestias Interact?</b> (Synestias Must Conserve Angular Momentum)

There are a few smaller, bite-sized notebooks that define terms or concepts that are commonly invoked in explanations throughout the main notebooks, namely:
* Angular Velocity   
* Keplerian Orbit 
* Midplane 
* Poincarè-Wavre Theorem
* Spherical and Cylindrical Coordinates 

A few potential Moon-forming synestia example cases (output from smoothed-particle hydrodynamics simulations) are used throughout the notebooks, especially as part of the interactive plots. Refer to `Synestia_Moon_Ex_Cases.ipynb` for details about these synestia example cases.

## Science Concepts Covered
These notebooks will touch on various science concepts: linear and angular velocity, linear and centripetal acceleration, conservation of angular momentum, gravitational forces acting on multiple bodies, non-conservative forces (gas drag), pressure gradients and hydrostatic equilibrium, material phases and phase changes, extreme thermal regimes for rocky materials, orbits (Keplerian), solid body rotation, giant impacts, deep time (early solar system), formation of moons and planets, and planetary structures (planet, moon, planetary disk, synestia).

## Student Learning Outcomes
In addition, these notebooks satisfy the following student learning outcomes based on CAS (Council for the Advancement of Standards in higher education) domains ([Wells, 2015](http://standards.cas.edu/getpdf.cfm?PDF=D87A29DC-D1D6-D014-83AA8667902C480B)).

1. Knowledge acquisition, construction, integration and application 
    * <i>Understanding knowledge from a range of disciplines</i>: "Possesses knowledge of...the physical world; possesses knowledge of [a specific] one or more subjects" (e.g. synestias, formation of Earth-Moon system)
    * <i>Constructing knowledge</i>: "Personalizes learning; makes meaning from text, instruction, and experience; uses experience and other sources of information to create new insights; ...recognizes one’s own capacity to create new understandings from learning activities"
    * <i>Relating knowledge to daily life</i>: "makes connections between classroom and out-of-classroom learning"


2. Cognitive complexity
    * <i>Critical thinking</i>: "analyzes, interprets, and makes judgments of the relevance and quality of information; assesses assumptions and considers alternative perspectives and solutions"  
    * <i>Reflective thinking</i>: "Applies previously understood information, concepts, and experiences to a new situation or setting; rethinks previous assumptions"  
    * <i>Effective reasoning</i>: "Uses complex information from a variety of sources including personal experience and observation to form a decision or opinion; is open to new ideas and perspectives"


3. Practical competence
    * <i>Technological competence</i>: "Demonstrates technological literacy and skills...stays current with technological innovations"

## Acknowledgements
I included previously published figures and videos throughout these notebooks. For each copyrighted material that appears in this thesis, I have referenced its respective source. The respective copyright owners have granted permission for their materials' use in this academic, non-profit publication. The materials have been reprinted with permission from: Donna Cox of the University of Illinois at Urbana-Champaign Visualization Center, The American Association for the Advancement of Science or AAAS (<i>Science</i>) and Robin M. Canup of Southwest Research Institute, American Geophysical Union or AGU (<i>Journal of Geophysical Research</i>) and Simon J. Lock, and Raluca Rufu. Email correspondence proving permission to use these copyrighted materials is located in the folder `copyright_permissions`.

## Contact
This set of Jupyter Notebooks is intended for non-profit educational purposes, in partial fulfillment of Gigja Hollyday's M.S. thesis in Geology at the University of California, Davis. For inquiries or feedback, please email gohollyday at ucdavis dot edu.

## References
Barba, L., & Forsyth, G. (2018). CFD python: The 12 steps to Navier-Stokes equations. <i>Journal of Open Source Education</i>, 1 (9), 21.

Freeman, S., Eddy, S. L., McDonough, M., Smith, M. K., Okoroafor, N., Jordt, H., & Wenderoth, M. P. (2014). Active learning increases student performance in science, engineering, and mathematics. <i>Proceedings of the National Academy of Sciences</i>, 111 (23), 8410-8415.

Kluyver, T., Ragan-Kelley, B., Pérez, F., Granger, B., Bussonnier, M., Frederic, J., Kelley, K., Hamrick, J., Grout, J., Corlay, S., Ivanov, P., Avila, D., Abdalla, S., and Willing, C. (2016). Jupyter notebooks -- a publishing format for reproducible computational work flows. In F. Loizides & B. Schmidt (Eds.), <i>Positioning and Power in Academic Publishing: Players, Agents and Agendas</i> (p. 87-90). IOS Press.

Lock, S. J., & Stewart, S. T. (2019). When Earth and the Moon were one. Retrieved from https://www.scientificamerican.com/article/when-earth-and-the-moon-were-one/ (Scientific American)

Rein, H., & Liu, S. F. (2012). REBOUND: an open-source multi-purpose N-body code for collisional dynamics. <i>Astronomy and Astrophysics</i>, 537 (A128), 1-10.

Rein, H., & Spiegel, D. S. (2015). IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits. <i>Monthly Notices of the Royal Astronomical Society</i>, 446 (2), 1424-1437.

Stewart, S. T. (2019). <i>Where did the Moon come from? A new theory.</i> TED Salon: U.S. Air Force. Retrieved from https://www.ted.com/talks/sarah_t_stewart_where_did_the_moon_come_from_a_new_theory?language=en (TED)

Wells, J. B. (Ed.). (2015). <i>CAS learning and development outcomes</i> (9th ed.). Council for the Advancement of Standards in Higher Education. Retrieved from http://standards.cas.edu/getpdf.cfm?PDF=D87A29DC-D1D6-D014-83AA8667902C480B (CAS)
