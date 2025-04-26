import { AtomicElementValue } from "./elements";
import { AtomicLabelValue } from "./labels";
export declare class ElementWithLabel {
    element: AtomicElementValue;
    label: AtomicLabelValue;
    constructor({ element, label }: {
        element: AtomicElementValue;
        label: AtomicLabelValue;
    });
    prettyPrint(): string;
}
