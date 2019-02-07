import assert from "assert";
import {Then} from "cucumber";

Then(/^Basis stored under "([^"]*)" key returns "([^"]*)" as formula$/, function (cacheKey, formula) {
    assert.equal(this["cacheKey"].formula, formula);
});
