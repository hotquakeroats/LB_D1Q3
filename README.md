# LB_D1Q3
Master's research code

Learning the physics of what I was doing was my primary goal.  As a result, please be advised that a great many standard software engineering practices went right out the window:

1. Globals are everywhere!  To readers with a strong computer science background, I'm sure you've gone veklempt.  I'll pause so you may collect yourselves... ... ... ... ... (younger readers, Youtube early-90s SNL sketches called "Coffee Talk" with Mike Myers as Linda Richmond)... ... ... ... better?  There is a time and a place for such things, and non-production code that is developed to better my physical understanding is both.  It's so much simpler to not have to trace a thousand setters and getters throughout a code base when you're learning such things, but please feel free to refactor as you see fit.
2. Security was an afterthought - buffer overflows, memory allocation/deallocation, fscanf vs. fgets, testing for NULL pointers, etc. were only addressed very superficially (if at all).  If you want to use this code as a base for anything more production-oriented, you'll really want to take a look at tightening the code up.
3. There are a number of places where a more computer sciencey algorithm would speed things up greatly (the elephant in the room being the brute force minimization).  I found that I developed a better intuition for the physics by eschewing such things, so I sacrificed clock cycles on purpose to get a better understanding.
4. Edge cases... I tended to only address the ones that specifically cause me problems.  There are more out there.
