# ddg-exercises-js

### Forked by Joseph Adrian
This fork will be updated as I complete the implementation of each assignment in the `projects` folder.
This is being done as part of an undergraduate independent study project at Towson University in Spring 2022.
Written below is a link to the course I am following as part of my independent study, as well as a link to the 
original repository for the skeleton code, followed by the original unaltered content of the README.md file from
the aformentioned repository. 
Credit goes to Keenan Crane as the instructor of the course, and Rohan Sawhney as the author of the original code.
[Course Link](https://brickisland.net/DDGSpring2022/)
[Original Repo](https://github.com/cmu-geometry/ddg-exercises-js)

ddg-exercises-js is a fast and flexible framework for 3D geometry processing
on the web! Easy integration with HTML/WebGL makes it particularly suitable for
things like mobile apps, online demos, and course content. For many tasks,
performance comes within striking distance of native (C++) code. Plus, since the
framework is pure JavaScript, **no compilation or installation** is necessary on any
platform. Moreover, geometry processing algorithms can be **edited in the browser**
(using for instance the [JavaScript Console](https://developers.google.com/web/tools/chrome-devtools/console/) in Chrome).

At a high level, the framework is divided into three parts - an implementation of
a halfedge mesh data structure, an optimized linear algebra package and skeleton
code for various geometry processing algorithms. Each algorithm comes with its own
viewer for rendering.

Detailed documentation and unit tests for each of these parts can be found in the docs
and tests directories of this [repository](https://github.com/cmu-geometry/ddg-exercises-js).

## Getting started

1. Clone the repository and change into the projects directory
```
git clone https://github.com/cmu-geometry/ddg-exercises-js.git
cd ddg-exercises-js/projects
```

2. Open the index.html file in any of the sub directories in a browser of your choice
(Chrome and Firefox usually provide better rendering performance than Safari).

## Dependencies (all included)

1. Linear Algebra - A wrapper around the C++ library [Eigen](https://eigen.tuxfamily.org) compiled
to [asm.js](http://asmjs.org) with [emscripten](http://emscripten.org). Future updates will compile
the more optimized sparse matrix library [Suitesparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) to asm.js.

2. Rendering - [three.js](https://threejs.org)

3. Unit Tests - [Mocha](http://mochajs.org) and [Chai](http://chaijs.com)

## About Javascript

The implementation of ddg-exercises-js attempts to minimize the use of obscure
Javascript language features. It should not be too difficult for anyone with experience
in a dynamic language like Python or familiar with the principles of Object Oriented Programming
to get a handle on Javascript syntax by reading through some of the code in this framework.
The documentation contains examples specific to this framework which will also be of help.
For a more formal introduction to Javascript, checkout this really nice [tutorial](https://javascript.info).

## Author

Rohan Sawhney

Email: rohansawhney@cs.cmu.edu

## License

[MIT](https://opensource.org/licenses/MIT)
