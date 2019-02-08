import {When} from "cucumber";

import {parseTable} from "../../table";
import {Basis} from "../../../../../src/basis/basis";

When(/^Basis is created with the following data:$/, function (table) {
    const config = parseTable(table, this)[0];
    this["cacheKey"] = new Basis(config.basis);
});
