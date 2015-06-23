---
layout: post
title: "When `True` is not `True`"
---

Today, whilst continuing development on [`Goldilocks`](https://github.com/SamStudio8/goldilocks), I
discovered a minor oddity that left me a little confused and bemused before lunch: `True` did not appear to be `True`...

Part of `Goldilocks`' functionality allows for the filtering of results; users may specify a dictionary of criteria
whose keys map to functions to be applied to result sets to perform such filtering. For example, `start_lte` and
`start_gte` both call the following function (with `operand` set -1 and 1, respectively), filtering out regions
whose starting base position is less than or greater than or equal to some value:

```python
def __exclude_start(region_dict, operand, position):
    if operand < 0:
        return region_dict["pos_start"] <= position
    elif operand > 0:
        return region_dict["pos_start"] >= position
    return False
```

Each of these exclusion checking functions are expected to return `True` if the criteria for exclusion has been met.
Following each different criteria check on the current region, the following mops up to see whether
the comparisons can be aborted early:

```python
[...]
elif name == "start_lte":
    ret = __exclude_start(region_dict, -1, to_apply["start_lte"]
[...]

if use_and:
    # Require all exclusions to be true... 
    if ret is False:
        return False
else:
    if ret is True:
        # If we're not waiting on all conditions, we can exclude on the first
        return True
```

However, during some testing this morning, I noticed spurious results: regions that I expected to be
excluded were not. The test suite confirms. I played around with the direction of my angle brackets
and switched around `True` and `False` to no avail.

I added a simple `print` statement to `__exclude_start`, everything appeared to behave as expected - `True` and `False`
were being printed as one would expect for each pair of positions. Yet the clean-up `if ret is True`
statement was definitely being "ignored".

```python
print("%d <= %d: %r" % (region_dict["pos_start"], position, region_dict["pos_start"] <= position))
>>> 1 <= 2: True
>>> 2 <= 1: False
```

Struggling for ideas I thought: maybe I'm not supposed to be using `return` like that?
I reduced the exclusion testing function and "manually" returned `True` where necessary:

```python
def __exclude_lte(a, b):
    if a <= b:
        return True
    return False
```

The test suite passes.

I try something else.

```python
def __exclude_lte(a, b):
    return bool(a <= b)
```

The test suite passes. What weird funky type magic is happening? I'm pretty certain I'm allowed to use `return`
like this and expressions should be automatically `bool` anyway?

```python
print("%d <= %d: %r (%s)" % (a, b, a <= b, type(a <= b)))
>>> 1 <= 2: True (<type 'numpy.bool_'>)
>>> 2 <= 1: False (<type 'numpy.bool_'>)
```

Oh crumbs. Now it all makes sense...

`Goldilocks` makes extensive use of the `numpy` package (primarily for its nice arrays)
which apparently implements its own boolean type that is returned when forming expressions that involve
other `numpy` types, such as `int64`. In Python, the `is` operator checks whether or not two variables
point at the same object in memory, it does not check for equality. Of course, here: `True` is not `numpy.bool_(True)`[^1]
and this is why `if ret is True` failed and results were not filtered.

Of course, as usual this is all my fault and could have been easily avoided. The anal C programmer in me likes explicit
checking of these sorts of things (and `is (not) None` is a frequent occurrence in my Python scripts)
but this whole trouble would have been avoided if I'd just used appropriate Python style and ditched the redundant
parts of the clean-up statements anyway:

```python
if use_and:
    # Require all exclusions to be true... 
    if not ret:
        return False
else:
    if ret:
        # If we're not waiting on all conditions, we can exclude on the first
        return True
```

* * *

# tl;dr
* TIL: `numpy` has its own `bool` type.
* One should be careful to remember the difference between testing identity (`is`) and equality (`==`) in Python.
* One should probably be more careful to avoid problems like this in the first place by using the language constructs properly...

[^1] Although:
    
    ```python
    True is not np.bool(True)
    >>> False
    ```
