# ARC

The following performance optimisations have been implemented:
* The usual store/restore of intermediate values.
* re-compute brances only if any of the inputs have changed.
* Recompute rates corresponding to branches where the corresponding sampling probability (categories) have changed.
* Recompute rates for branches with length changes.

Issues:
* There is an over-reporting of changes to the omega parameter when using param.somethingIsDirty(). Currently investigating.

