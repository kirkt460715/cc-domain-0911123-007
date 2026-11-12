"""Exports for method agreement studies."""

import csv
import io

from ..values import canonical


class AgreementExports:
    def __init__(self, studies):
        self.studies = studies

    def json_bytes(self, study_id):
        study = self.studies.get(study_id)
        return (canonical(study["payload"]) + "\n").encode("utf-8")

    def pair_rows(self, study_id):
        """One CSV row per paired reading (replicates stay separate)."""
        payload = self.studies.get(study_id)["payload"]
        output = io.StringIO(newline="")
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(
            [
                "specimen",
                "replicate",
                "reference",
                "candidate",
                "difference",
                "average",
                "reference_observation",
                "candidate_observation",
                "reference_revision",
                "candidate_revision",
                "reference_certificate",
                "candidate_certificate",
                "reference_at",
                "candidate_at",
            ]
        )
        for point in payload["analysis"].get("points", []):
            writer.writerow(
                [
                    point["specimen"],
                    point["replicate"],
                    point["reference"],
                    point["candidate"],
                    point["difference"],
                    point["average"],
                    point["reference_observation"],
                    point["candidate_observation"],
                    point["reference_revision"],
                    point["candidate_revision"],
                    point["reference_certificate"],
                    point["candidate_certificate"],
                    point["reference_at"],
                    point["candidate_at"],
                ]
            )
        return output.getvalue()

    def exclusion_rows(self, study_id):
        payload = self.studies.get(study_id)["payload"]
        output = io.StringIO(newline="")
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(["code", "side", "observation", "message", "evidence"])
        for exclusion in payload["pairing"]["exclusions"]:
            writer.writerow(
                [
                    exclusion["code"],
                    exclusion["side"],
                    exclusion["observation"],
                    exclusion["message"],
                    canonical(exclusion["evidence"]),
                ]
            )
        return output.getvalue()

    def summary_rows(self, study_id):
        payload = self.studies.get(study_id)["payload"]
        analysis = payload["analysis"]
        output = io.StringIO(newline="")
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(["metric", "value"])
        ba = analysis.get("bland_altman")
        if ba:
            for key in (
                "pairs",
                "specimens",
                "scale",
                "bias",
                "sd_difference",
                "lower_limit",
                "upper_limit",
                "correlation",
                "repeatability_limit",
            ):
                writer.writerow([key, ba[key]])
        decision = analysis["decision"]
        writer.writerow(["decision", decision["status"]])
        writer.writerow(["reasons", ";".join(decision["reasons"])])
        writer.writerow(
            ["acceptable_difference", decision["acceptable_difference"]]
        )
        if analysis.get("deming"):
            deming = analysis["deming"]
            writer.writerow(["deming_slope", deming["slope"]])
            writer.writerow(["deming_intercept", deming["intercept"]])
            writer.writerow(["deming_variance_ratio", deming["variance_ratio"]])
        return output.getvalue()
