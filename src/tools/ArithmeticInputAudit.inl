// This input check deliberately stops before constructing binary64 solver state.
std::string arithmeticInputHex(vela::ArithmeticReal value)
{
    if (!std::isfinite(value)) throw std::invalid_argument("Arithmetic input conversion overflowed");
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << std::hexfloat << value;
    return stream.str();
}

nlohmann::json checkMaterialArithmetic(const nlohmann::json& material, vela::ArithmeticReal temperature,
                                        bool silicon)
{
    const auto values = vela::evaluateSemiconductorArithmetic(material, temperature);
    nlohmann::json result{{"gap_eV_hex", arithmeticInputHex(values.bandgap_eV)},
                          {"nc_cm3_hex", arithmeticInputHex(values.electronsDOS_cm3)},
                          {"nv_cm3_hex", arithmeticInputHex(values.holesDOS_cm3)},
                          {"ni_cm3_hex", arithmeticInputHex(values.intrinsicDensity_cm3)}};
    if (silicon)
        result["reference_phi_V_hex"] = arithmeticInputHex(values.referencePotential_V);
    return result;
}

nlohmann::json checkArithmeticInput(const std::string& path)
{
    std::ifstream file(path);
    if (!file) throw std::invalid_argument("Cannot open arithmetic input file");
    nlohmann::json input;
    file >> input;
    if (input.at("schema") != "vela.arithmetic_input.v1")
        throw std::invalid_argument("Unsupported arithmetic input schema");
    vela::constants::requireProfile(input.at("physical_constants_profile").get<std::string>());
    const auto temperature = vela::parseArithmeticDecimal(input.at("temperature_K").get<std::string>());
    if (temperature <= 0) throw std::invalid_argument("Arithmetic input temperature must be positive");
    auto output = arithmeticCapabilities();
    output["simulation_started"] = false;
    output["temperature_K_hex"] = arithmeticInputHex(temperature);
    output["thermal_voltage_V_hex"] = arithmeticInputHex(vela::constants::thermalVoltage(temperature));
    if (input.contains("materials")) {
        for (const auto& [name, material] : input.at("materials").items())
            output["materials"][name] = checkMaterialArithmetic(material, temperature, name == "Si");
    }
    for (const auto& [role, source] : input.at("sources").items()) {
        auto& check = output["sources"][role];
        check["sha256"] = source.at("sha256");
        const bool binary64 = source.value("number_encoding", std::string{"decimal"}) == "binary64_roundtrip";
        if (binary64) {
            check["number_encoding"] = "binary64_roundtrip";
            check["original_tdr_sha256"] = source.at("original_tdr_sha256");
        }
        auto& values = check["values_hex"];
        values = nlohmann::json::object();
        for (const auto& [key, literal] : source.at("numbers").items()) {
            try {
                const auto text = literal.get<std::string>();
                const auto value = binary64
                    ? static_cast<vela::ArithmeticReal>(vela::parseArithmeticDecimal<double>(text))
                    : vela::parseArithmeticDecimal(text);
                values[key] = arithmeticInputHex(value);
                if (role == "grd" || role == "tdr_coordinates")
                    check["meters_hex"][key] = arithmeticInputHex(value * VELA_ARITHMETIC_LITERAL(1e-6));
                else if (role == "dat" || role == "tdr_doping")
                    check["per_m3_hex"][key] = arithmeticInputHex(value * VELA_ARITHMETIC_LITERAL(1e6));
            } catch (const std::exception& error) {
                throw std::invalid_argument(role + " entry " + key + ": " + error.what());
            }
        }
        check["number_count"] = values.size();
    }
    return output;
}

nlohmann::json verifyExtended80StageState(std::istream& request)
{
    if constexpr (vela::compiledArithmeticPrecision != vela::ArithmeticPrecision::Extended80)
        throw std::invalid_argument("Stage state requires an extended80 runner");
    const auto contract = nlohmann::json::parse(request);
    if (contract.at("schema") != "vela.verify-extended80-stage-state.v1")
        throw std::invalid_argument("Unsupported extended80 stage-state request schema");
    const auto& nodes = contract.at("node_count");
    if (!nodes.is_number_unsigned() || nodes.get<std::uint64_t>() == 0)
        throw std::invalid_argument("Stage state requires a nonempty mesh");
    const vela::ExtendedPopulationArchiveIdentity identity{
        contract.at("mesh_sha256").get<std::string>(),
        contract.at("input_sha256").get<std::string>(),
        vela::ExtendedPopulationArchiveRole::Accepted};
    const auto state = vela::readExtendedPopulationArchive(
        contract.at("archive").get<std::string>(), nodes.get<std::uint64_t>(), identity);
    const auto exact = [&](const char* key) {
        return vela::parseArithmeticDecimal<long double>(contract.at(key).get<std::string>());
    };
    if (state.deviceBias_V != exact("device_voltage_V") ||
        state.supplyBias_V != exact("supply_voltage_V") ||
        state.drainCurrent_A_per_um != exact("current_A_per_um"))
        throw std::invalid_argument("Stage voltage or current differs from the authoritative archive");
    return {{"verified", true}, {"arithmetic_precision", "extended80"}};
}
